"""Loads everything the replica needs from the vendored NASA sources.

* MWS workspace from SimSetup/setup_*.m (real maps, via the mini MATLAB evaluator)
* mask parameters and block wiring from AGTF30_eng.mdl (text Simulink model)
* gas-table coefficients and constants from the T-MATS C sources
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import numpy as np

from ._matlab import Evaluator, MatlabError, read_text, run_script, to_python


# ------------------------------------------------------------------ .mdl text parser
def parse_mdl(path: str | Path) -> dict:
    """Parse a text .mdl file into nested nodes {'t': type, 'p': params, 'k': children}."""
    root = {"t": "root", "p": {}, "k": []}
    stack = [root]
    last = None
    for ln in read_text(path).split("\n"):
        s = ln.strip()
        if not s:
            continue
        m = re.match(r"^([A-Za-z_$][\w.$]*)\s*\{$", s)
        if m:
            node = {"t": m.group(1), "p": {}, "k": []}
            stack[-1]["k"].append(node)
            stack.append(node)
            last = None
        elif s == "}":
            stack.pop()
            last = None
        elif s.startswith('"') and last is not None:
            stack[-1]["p"][last] += s[1:-1] if s.endswith('"') else s[1:]
        else:
            m = re.match(r"^(\S+)\s+(.*)$", s)
            if m:
                v = m.group(2)
                if v.startswith('"'):
                    v = v[1:-1] if v.endswith('"') and len(v) > 1 else v[1:]
                stack[-1]["p"][m.group(1)] = v
                last = m.group(1)
    return root


def _kids(node: dict, typ: str) -> list[dict]:
    return [k for k in node["k"] if k["t"] == typ]


def _name(block: dict) -> str:
    return block["p"].get("Name", "").replace("\\n", " ")


class EngineModel:
    """Block parameters and resolved wiring of AGTF30_eng.mdl."""

    def __init__(self, path: str | Path):
        root = parse_mdl(path)
        model = _kids(root, "Model")[0]
        self.defaults: dict[str, dict] = {}
        for d in _kids(model, "BlockParameterDefaults"):
            for b in _kids(d, "Block"):
                self.defaults[b["p"].get("BlockType", "")] = b["p"]
        self.blocks: dict[str, dict] = {}
        self.systems: dict[str, dict] = {}
        self._walk(_kids(model, "System")[0], "")

    def _walk(self, system: dict, path: str) -> None:
        self.systems[path] = system
        for b in _kids(system, "Block"):
            p = (path + "/" if path else "") + _name(b)
            params = dict(b["p"])
            for obj in _kids(b, "Object"):                 # mask parameters of a masked subsystem
                if obj["p"].get("$ClassName") == "Simulink.Mask":
                    for arr in _kids(obj, "Array"):
                        for mp in _kids(arr, "Object"):
                            if "Name" in mp["p"]:
                                params[mp["p"]["Name"]] = mp["p"].get("Value", "")
            self.blocks[p] = params
            for s in _kids(b, "System"):
                self._walk(s, p)

    def param(self, block: str, key: str) -> str:
        """Raw string of a block (mask) parameter, falling back to the block-type default."""
        b = self.blocks[block]
        if key in b:
            return b[key]
        return self.defaults.get(b.get("BlockType", ""), {})[key]

    def connections(self, system: str) -> set[tuple[str, int, str, int]]:
        """All (src block, src port, dst block, dst port) of a system, Goto/From tags resolved."""
        sysn = self.systems[system]
        raw: list[tuple[str, int, str, int]] = []

        def rec(node: dict, src: tuple[str, int]) -> None:
            p = node["p"]
            if "DstBlock" in p:
                port = p["DstPort"]
                raw.append((src[0], src[1], p["DstBlock"].replace("\\n", " "),
                            int(port) if port.isdigit() else 0))
            for br in _kids(node, "Branch"):
                rec(br, src)

        for ln in _kids(sysn, "Line"):
            rec(ln, (ln["p"]["SrcBlock"].replace("\\n", " "), int(ln["p"]["SrcPort"])))
        btype, tag = {}, {}
        for b in _kids(sysn, "Block"):
            btype[_name(b)] = b["p"].get("BlockType")
            tag[_name(b)] = b["p"].get("GotoTag")
        feeds = {tag[d]: (s, sp) for (s, sp, d, dp) in raw if btype.get(d) == "Goto"}
        out = set()
        for (s, sp, d, dp) in raw:
            if btype.get(d) == "Goto":
                continue
            if btype.get(s) == "From":
                if tag[s] not in feeds:
                    continue
                s, sp = feeds[tag[s]]
            out.add((s, sp, d, dp))
        return out

    def port_of(self, system: str, block: str) -> int:
        """Port number of an Inport/Outport block inside a subsystem."""
        return int(self.blocks[system + "/" + block].get("Port", "1"))


# Wiring the replica's engine_replica.Engine.residuals() assumes (checked at load time).
EXPECTED_WIRING: dict[str, list[tuple[str, int, str, int]]] = {
    "Engine": [
        ("Ind", 1, "Demux", 1), ("Demux", 1, "Ambient", 1), ("Demux", 2, "Core", 4),
        ("Demux", 3, "Core", 3), ("Demux", 4, "FrontEng", 3), ("Demux", 5, "FrontEng", 4),
        ("Demux", 6, "FrontEng", 5), ("Demux", 7, "BackEng", 4),
        ("FrontEng", 3, "Mux2", 1), ("FrontEng", 4, "Mux2", 2), ("Core", 3, "Mux2", 3),
        ("Core", 4, "Mux2", 4), ("BackEng", 1, "Mux2", 5), ("BackEng", 2, "Mux2", 6),
        ("Bypass", 1, "Mux2", 7), ("Mux2", 1, "Dep", 1),
        ("Sys_In", 1, "Demux2", 1), ("Demux2", 1, "Ambient", 2), ("Demux2", 2, "Ambient", 3),
        ("Demux2", 3, "Ambient", 4), ("Demux2", 4, "Core", 1), ("Demux2", 5, "Bypass", 1),
        ("Demux2", 6, "FrontEng", 7),
        ("Nmech", 1, "Demux1", 1), ("Demux1", 1, "LP_Shaft", 3), ("Demux1", 2, "HP_Shaft", 3),
        ("LP_Shaft", 1, "FrontEng", 6), ("LP_Shaft", 1, "BackEng", 5), ("HP_Shaft", 1, "Core", 5),
        ("FrontEng", 6, "Mux", 1), ("FrontEng", 5, "Mux", 2), ("BackEng", 3, "Product", 1),
        ("Constant", 1, "Product", 2), ("Product", 1, "Mux", 3), ("Mux", 1, "LP_Shaft", 1),
        ("Constant1", 1, "LP_Shaft", 2),
        ("Core", 5, "Mux1", 1), ("Core", 6, "Mux1", 2), ("Mux1", 1, "HP_Shaft", 1),
        ("PwrDebit [HP]", 1, "HP_Shaft", 2),
        ("Core", 1, "BackEng", 1), ("Core", 2, "BackEng", 2), ("FrontEng", 1, "Bypass", 3),
        ("FrontEng", 2, "Core", 2), ("Ambient", 1, "FrontEng", 1), ("Ambient", 2, "FrontEng", 2),
        ("Ambient", 2, "Bypass", 2), ("Ambient", 2, "BackEng", 3),
        ("LP_Shaft", 2, "Mux3", 1), ("HP_Shaft", 2, "Mux3", 2), ("Mux3", 1, "Ndot", 1),
    ],
    "Engine/FrontEng": [
        ("AmbFlow", 1, "Inlet", 1), ("Pamb", 1, "Inlet", 2), ("Inlet", 1, "FAN", 1),
        ("Fan_Rline", 1, "FAN", 2), ("GearBox", 1, "FAN", 3), ("FAN", 5, "GearBox", 1),
        ("LP_Nmech", 1, "GearBox", 2), ("GearBox", 2, "Fan_Trq", 1), ("FAN", 3, "Splitter", 1),
        ("Splitter_BPR", 1, "Splitter", 2), ("Splitter", 2, "duct2", 1), ("duct2", 1, "LPC", 1),
        ("LPC_Rline", 1, "LPC", 2), ("LP_Nmech", 1, "LPC", 3), ("LPC", 3, "Valve", 3),
        ("Splitter", 1, "Valve", 1), ("VBVdmd", 1, "Valve", 2), ("Valve", 2, "duct25", 1),
        ("Valve", 1, "BypFlow", 1), ("duct25", 1, "FlwOut", 1), ("FAN", 4, "NErr1", 1),
        ("LPC", 4, "NErr2", 1), ("LPC", 5, "LPC_Trq", 1),
    ],
    "Engine/Core": [
        ("FlwIn", 1, "HPC", 1), ("HPC_Rline", 1, "HPC", 2), ("HP_Nmech", 1, "HPC", 3),
        ("HPC", 3, "Burner", 2), ("Wf", 1, "Burner", 1), ("Burner", 1, "HPT", 2),
        ("HPC", 2, "Demux", 1), ("Demux", 1, "Bld", 1), ("Demux", 2, "Mux", 1),
        ("Demux", 3, "Mux", 2), ("Mux", 1, "HPT", 1), ("HPT_PRIn", 1, "HPT", 3),
        ("HP_Nmech", 1, "HPT", 4), ("HPT", 1, "FlwOut", 1), ("HPC", 4, "NErr3", 1),
        ("HPT", 2, "NErr4", 1), ("HPC", 5, "HPC_Trq", 1), ("HPT", 3, "HPT_Trq", 1),
        ("HPC", 3, "S36 Static Calc", 1),
    ],
    "Engine/BackEng": [
        ("FlowIn", 1, "duct45", 1), ("duct45", 1, "LPT", 2), ("BldFlw", 1, "LPT", 1),
        ("LPT_PRIn", 1, "LPT", 3), ("LP_Nmech", 1, "LPT", 4), ("LPT", 1, "duct5", 1),
        ("duct5", 1, "CoreNoz", 1), ("Pamb", 1, "CoreNoz", 2), ("LPT", 2, "NErr5", 1),
        ("CoreNoz", 3, "NErr6", 1), ("LPT", 3, "LPT_Trq", 1),
    ],
    "Engine/Bypass": [
        ("BypFlow", 1, "duct17", 1), ("duct17", 1, "BypNoz", 1), ("Pamb", 1, "BypNoz", 2),
        ("VAFNdmd", 1, "BypNoz", 3), ("VAFNdmd", 1, "BypNoz", 4), ("BypNoz", 3, "NErr7", 1),
    ],
}
EXPECTED_PORTS: dict[str, dict[str, int]] = {
    "Engine/FrontEng": {"AmbFlow": 1, "Pamb": 2, "Fan_Rline": 3, "Splitter_BPR": 4, "LPC_Rline": 5,
                        "LP_Nmech": 6, "VBVdmd": 7, "BypFlow": 1, "FlwOut": 2, "NErr1": 3,
                        "NErr2": 4, "Fan_Trq": 5, "LPC_Trq": 6},
    "Engine/Core": {"Wf": 1, "FlwIn": 2, "HPC_Rline": 3, "HPT_PRIn": 4, "HP_Nmech": 5, "Bld": 1,
                    "FlwOut": 2, "NErr3": 3, "NErr4": 4, "HPC_Trq": 5, "HPT_Trq": 6},
    "Engine/BackEng": {"BldFlw": 1, "FlowIn": 2, "Pamb": 3, "LPT_PRIn": 4, "LP_Nmech": 5,
                       "NErr5": 1, "NErr6": 2, "LPT_Trq": 3},
    "Engine/Bypass": {"VAFNdmd": 1, "Pamb": 2, "BypFlow": 3, "NErr7": 1},
    "Engine": {"Sys_In": 1, "Ind": 2, "Nmech": 3, "Eng_Out": 1, "Ndot": 2, "Dep": 3},
}


def check_wiring(model: EngineModel) -> list[str]:
    """Return the list of wiring assumptions NOT found in the .mdl (empty = all confirmed)."""
    missing = []
    for system, exp in EXPECTED_WIRING.items():
        have = model.connections(system)
        for c in exp:
            if c not in have:
                missing.append(f"{system}: {c[0]}:{c[1]} -> {c[2]}:{c[3]}")
    for system, ports in EXPECTED_PORTS.items():
        for blk, num in ports.items():
            if model.port_of(system, blk) != num:
                missing.append(f"{system}: port number of {blk} is not {num}")
    return missing


# ------------------------------------------------------------------ C sources
def parse_c_tables(mex_dir: Path) -> dict[str, list[float]]:
    """Gas-table arrays embedded in pt2sc_TMATS.c / t2hc_TMATS.c (ITAB, TTTAB, A/B/C/D...)."""
    tabs: dict[str, list[float]] = {}
    for fn in ("t2hc_TMATS.c", "h2tc_TMATS.c", "pt2sc_TMATS.c", "sp2tc_TMATS.c"):
        text = re.sub(r"/\*.*?\*/", "", read_text(mex_dir / fn), flags=re.S)
        for m in re.finditer(r"(?:int|double)\s+(\w+)\[(\d+)\]\s*=\s*\{([^}]*)\}", text):
            vals = [float(x) for x in m.group(3).replace("\n", " ").split(",") if x.strip()]
            if len(vals) != int(m.group(2)):
                raise ValueError(f"{fn}: table {m.group(1)} has wrong length")
            if m.group(1) in tabs and tabs[m.group(1)] != vals:
                raise ValueError(f"table {m.group(1)} differs between C files ({fn})")
            tabs[m.group(1)] = vals
    return tabs


def parse_c_constants(mex_dir: Path) -> dict[str, float]:
    """#define constants of constants_TMATS.h."""
    out = {}
    for m in re.finditer(r"#define\s+(\w+)\s+([-+0-9.eE]+)", read_text(mex_dir / "constants_TMATS.h")):
        out[m.group(1)] = float(m.group(2))
    return out


def library_constants(lib_dir: Path) -> dict[str, float]:
    """Constants of pure-Simulink library blocks, read from the library files when possible."""
    out = {"shaft_gain": -5252.11, "ambient_g": 32.1740, "secs2hrs": 3600.0}   # values seen in the source
    try:
        txt = read_text(lib_dir / "Lib_Turbo_Shaft_TMATS.mdl")
        m = re.search(r'Name\s+"C_HP_PER_RPMtoFT_LBF".*?Gain\s+"([-0-9.eE]+)"', txt, re.S)
        if m:
            out["shaft_gain"] = float(m.group(1))
    except OSError:
        pass
    try:
        import zipfile
        with zipfile.ZipFile(lib_dir / "Lib_Turbo_Ambient_TMATS.slx") as z:
            xml = z.read("simulink/blockdiagram.xml").decode("utf-8", "replace")
        m = re.search(r'Name="g" SID="\d+">.*?<P Name="Value">([-0-9.eE]+)</P>', xml, re.S)
        if m:
            out["ambient_g"] = float(m.group(1))
    except (OSError, KeyError):
        pass
    return out


# ------------------------------------------------------------------ everything together
class Sources:
    """All parsed NASA data: MWS (python values), engine .mdl model, C tables/constants."""

    def __init__(self, vendor_dir: str | Path):
        v = Path(vendor_dir)
        agtf = v / "agtf30" / "AGTF30"
        ws: dict[str, Any] = {}
        run_script(agtf / "SimSetup" / "setup_AllEng.m", ws, strict=True)
        self.ctrl_skipped = run_script(agtf / "SimSetup" / "setup_Controller.m", ws, strict=False)
        self._ws = ws
        self.MWS = to_python(ws["MWS"])
        self.model = EngineModel(agtf / "AGTF30_eng.mdl")
        self.wiring_missing = check_wiring(self.model)
        mex = v / "tmats" / "TMATS_Library" / "MEX"
        self.tables = parse_c_tables(mex)
        self.const = parse_c_constants(mex)
        self.lib = library_constants(v / "tmats" / "TMATS_Library")

    def mask(self, block: str, key: str) -> Any:
        """Evaluate a mask parameter of a block of AGTF30_eng.mdl in the MWS workspace."""
        raw = self.model.param(block, key).strip()
        if raw in ("on", "off"):
            return 1.0 if raw == "on" else 0.0
        if raw == "Convergent":
            return 1.0
        if raw == "Convergent-Divergent":
            return 2.0
        try:
            return to_python(Evaluator({"MWS": self._ws["MWS"]}).eval(raw.replace("\\n", " ")))
        except MatlabError as exc:
            raise MatlabError(f"mask {block}:{key} = {raw!r}: {exc}") from exc

    def vec(self, block: str, key: str) -> list[float]:
        """Mask parameter as a flat list of floats."""
        return [float(x) for x in np.atleast_1d(np.asarray(self.mask(block, key), dtype=float)).ravel()]
