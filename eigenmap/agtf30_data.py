"""Read NASA's AGTF30 parameter files (MATLAB text) and T-MATS gas tables (C text).

The NASA repository stores the component maps and constants as MATLAB scripts
(`vendor/agtf30/AGTF30/SimSetup/setup_*.m`, old-Mac CR line endings) and the gas
property polynomials as C arrays (`vendor/tmats/TMATS_Library/MEX/*.c`).
Nothing is transcribed by hand: this module contains a very small interpreter for
the subset of MATLAB used in those files (matrix literals, ranges, transposes,
indexing with `end`, row/column assignment, a few functions) and a regular
expression reader for the C arrays.

Course topic: Engine model (not a course topic) - data plumbing
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np

VENDOR = Path(__file__).resolve().parent.parent / "vendor"
DATA_DIR = Path(__file__).resolve().parent / "data"

SETUP_FILES = ("setup_Controller.m", "setup_FAN.m", "setup_LPC.m", "setup_HPC.m",
               "setup_HPT.m", "setup_LPT.m", "setup_Duct.m", "setup_AllEng.m")


# --------------------------------------------------------------------------
# 1. A tiny MATLAB-subset interpreter
# --------------------------------------------------------------------------
_TOKEN = re.compile(r"""
    (?P<num>(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?)
  | (?P<name>[A-Za-z_]\w*(\.[A-Za-z_]\w*)*)
  | (?P<op>\.\*|\./|\.\^|[-+*/^'()\[\],;:])
  | (?P<ws>[ \t]+)
  | (?P<nl>\n)
""", re.X)


class _Tok:
    __slots__ = ("kind", "text", "space_before", "space_after")

    def __init__(self, kind: str, text: str, space_before: bool):
        self.kind, self.text, self.space_before, self.space_after = kind, text, space_before, False


def _tokenize(src: str) -> list[_Tok]:
    """Split MATLAB source into tokens, remembering the blanks around each token.

    (Blanks matter inside [ ]: `[a -b]` is two elements, `[a - b]` is one.)
    """
    toks: list[_Tok] = []
    pos, space = 0, False
    while pos < len(src):
        m = _TOKEN.match(src, pos)
        if not m:
            raise ValueError("MATLAB tokenizer: cannot read %r" % src[pos:pos + 30])
        pos = m.end()
        kind = m.lastgroup
        if kind == "ws":
            space = True
            if toks:
                toks[-1].space_after = True
            continue
        if kind == "nl":
            kind, text = "op", ";"          # a newline inside [ ] ends a row
            toks.append(_Tok(kind, "\n", space))
            space = True
            continue
        toks.append(_Tok(kind, m.group(kind), space))
        space = False
    return toks


class _Parser:
    """Recursive-descent evaluator for MATLAB expressions (values are floats or 2-D numpy arrays)."""

    def __init__(self, toks: list[_Tok], env: dict):
        self.t, self.i, self.env = toks, 0, env
        self.in_matrix = 0
        self.end_stack: list[int] = []

    # -- helpers
    def peek(self) -> _Tok | None:
        return self.t[self.i] if self.i < len(self.t) else None

    def take(self) -> _Tok:
        tok = self.t[self.i]
        self.i += 1
        return tok

    def accept(self, text: str) -> bool:
        tok = self.peek()
        if tok is not None and tok.kind == "op" and tok.text == text:
            self.i += 1
            return True
        return False

    def skip_newlines(self) -> None:
        while self.peek() is not None and self.peek().text == "\n":
            self.i += 1

    # -- grammar
    def expr(self):
        first = self.add()
        tok = self.peek()
        if tok is not None and tok.text == ":" and not self._colon_is_index():
            self.take()
            second = self.add()
            if self.accept(":"):
                third = self.add()
                step, stop = float(second), float(third)
            else:
                step, stop = 1.0, float(second)
            n = int(np.floor((stop - float(first)) / step + 1e-10)) + 1
            return (float(first) + step * np.arange(n)).reshape(1, -1)
        return first

    def _colon_is_index(self) -> bool:
        return False

    def _binary_minus_allowed(self, tok: _Tok) -> bool:
        """Inside [ ], `a -b` starts a new element, `a - b` and `a-b` are a subtraction."""
        if not self.in_matrix:
            return True
        return not (tok.space_before and not tok.space_after)

    def add(self):
        val = self.mul()
        while True:
            tok = self.peek()
            if tok is None or tok.kind != "op" or tok.text not in "+-" or not self._binary_minus_allowed(tok):
                return val
            self.take()
            rhs = self.mul()
            val = val + rhs if tok.text == "+" else val - rhs

    def mul(self):
        val = self.unary()
        while True:
            tok = self.peek()
            if tok is None or tok.kind != "op" or tok.text not in ("*", "/", ".*", "./"):
                return val
            self.take()
            rhs = self.unary()
            if tok.text in ("*", ".*"):
                if tok.text == "*" and isinstance(val, np.ndarray) and isinstance(rhs, np.ndarray):
                    val = val @ rhs
                else:
                    val = val * rhs
            else:
                val = val / rhs

    def unary(self):
        tok = self.peek()
        if tok is not None and tok.kind == "op" and tok.text in "+-":
            self.take()
            val = self.unary()
            return -val if tok.text == "-" else val
        return self.power()

    def power(self):
        base = self.postfix()
        tok = self.peek()
        if tok is not None and tok.kind == "op" and tok.text in ("^", ".^"):
            self.take()
            return base ** self.unary()
        return base

    def postfix(self):
        val = self.primary()
        while True:
            tok = self.peek()
            if tok is None or tok.kind != "op":
                return val
            if tok.text == "'" and not tok.space_before:
                self.take()
                val = val.T if isinstance(val, np.ndarray) else val
            elif tok.text == "(" and not tok.space_before and isinstance(val, np.ndarray):
                self.take()
                val = self.index(val)
            else:
                return val

    def index(self, arr: np.ndarray):
        """arr(i), arr(i,j), arr(:,j), arr(end), arr(end - 1) with MATLAB 1-based indices."""
        saved, self.in_matrix = self.in_matrix, 0
        args = []
        dim = 0
        while True:
            tok = self.peek()
            if tok.text == ":" and self.t[self.i + 1].text in (",", ")"):
                self.take()
                args.append(slice(None))
            else:
                nargs_total = 1 if self._single_index_ahead() else 2
                self.end_stack.append(arr.size if nargs_total == 1 else arr.shape[dim])
                v = self.expr()
                self.end_stack.pop()
                args.append(int(round(float(v))) - 1)
            dim += 1
            if self.accept(","):
                continue
            if self.accept(")"):
                break
            raise ValueError("MATLAB index: expected , or )")
        self.in_matrix = saved
        if len(args) == 1:
            flat = arr.reshape(-1, order="F")
            out = flat[args[0]]
        else:
            out = arr[args[0], args[1]]
        if isinstance(out, np.ndarray):
            if out.ndim == 1:      # MATLAB keeps orientation: X(:,j) is a column, X(i,:) a row
                out = out.reshape(-1, 1) if isinstance(args[0], slice) else out.reshape(1, -1)
            return out
        return float(out)

    def _single_index_ahead(self) -> bool:
        depth, j = 0, self.i
        while j < len(self.t):
            x = self.t[j].text
            if x in "([":
                depth += 1
            elif x in ")]":
                if depth == 0:
                    return True
                depth -= 1
            elif x == "," and depth == 0:
                return False
            j += 1
        return True

    def primary(self):
        tok = self.take()
        if tok.kind == "num":
            return float(tok.text)
        if tok.kind == "name":
            if tok.text == "end":
                return float(self.end_stack[-1])
            if tok.text == "pi":
                return float(np.pi)
            nxt = self.peek()
            if tok.text in _FUNCS and nxt is not None and nxt.text == "(":
                self.take()
                saved, self.in_matrix = self.in_matrix, 0
                args = []
                if not self.accept(")"):
                    while True:
                        args.append(self.expr())
                        if self.accept(")"):
                            break
                        if not self.accept(","):
                            raise ValueError("MATLAB call: expected , or )")
                self.in_matrix = saved
                return _FUNCS[tok.text](*args)
            if tok.text not in self.env:
                raise KeyError("MATLAB name not defined: %s" % tok.text)
            return self.env[tok.text]
        if tok.text == "(":
            saved, self.in_matrix = self.in_matrix, 0
            val = self.expr()
            self.in_matrix = saved
            if not self.accept(")"):
                raise ValueError("MATLAB: expected )")
            return val
        if tok.text == "[":
            return self.matrix()
        raise ValueError("MATLAB: unexpected token %r" % tok.text)

    def matrix(self):
        """[a b; c d] : blanks/commas separate elements, ; or newline separate rows."""
        self.in_matrix += 1
        rows, row = [], []
        while True:
            tok = self.peek()
            if tok is None:
                raise ValueError("MATLAB: unterminated [")
            if tok.text == "]":
                self.take()
                break
            if tok.text in (";", "\n"):
                self.take()
                if row:
                    rows.append(row)
                row = []
                continue
            if tok.text == ",":
                self.take()
                continue
            row.append(self.expr())
        if row:
            rows.append(row)
        self.in_matrix -= 1
        if not rows:
            return np.zeros((0, 0))
        blocks = [np.hstack([np.atleast_2d(np.asarray(e, dtype=float)) for e in r]) for r in rows]
        return np.vstack(blocks)


def _ones(*shape):
    return np.ones(tuple(int(s) for s in shape)) if len(shape) > 1 else np.ones((int(shape[0]), int(shape[0])))


def _max(a, b=None):
    if b is not None:
        return np.maximum(a, b)
    a = np.atleast_2d(a)
    return float(np.max(a)) if 1 in a.shape else np.max(a, axis=0, keepdims=True)


def _min(a, b=None):
    if b is not None:
        return np.minimum(a, b)
    a = np.atleast_2d(a)
    return float(np.min(a)) if 1 in a.shape else np.min(a, axis=0, keepdims=True)


_FUNCS = {
    "ones": _ones,
    "zeros": lambda *s: np.zeros(tuple(int(v) for v in s)),
    "max": _max, "min": _min,
    "length": lambda a: float(max(np.atleast_2d(a).shape)),
    "numel": lambda a: float(np.atleast_2d(a).size),
    "sqrt": np.sqrt,
}

_LHS = re.compile(r"^\s*([A-Za-z_]\w*(?:\.[A-Za-z_]\w*)*)\s*(\(([^=]*)\))?\s*=(?!=)")


def _statements(text: str) -> list[str]:
    """Normalise line endings, strip comments, join `...` continuations, split statements."""
    text = text.replace("\r\n", "\n").replace("\r", "\n")
    lines = []
    for line in text.split("\n"):
        # strip comments: a % that is not inside a quoted string (these files have no '%' strings)
        cut = line.find("%")
        if cut >= 0:
            line = line[:cut]
        lines.append(line.rstrip())
    joined, buf = [], ""
    for line in lines:
        if line.endswith("..."):
            buf += line[:-3] + " "
            continue
        buf += line
        joined.append(buf)
        buf = ""
    # a statement may span several lines when a [ ... ] is open (rows separated by newlines)
    out, cur, depth = [], "", 0
    for line in joined:
        cur = (cur + "\n" + line) if cur else line
        depth += line.count("[") - line.count("]")
        if depth <= 0:
            if cur.strip():
                out.append(cur.strip())
            cur, depth = "", 0
    return out


def run_matlab(text: str, env: dict, skipped: list[str]) -> None:
    """Execute the assignments of one setup file into `env` (flat dict: 'MWS.FAN.NcVec' -> value)."""
    for st in _statements(text):
        st = st.rstrip(";").strip()
        if not st or st.startswith(("function", "end", "switch", "case", "otherwise", "clear", "if ", "else")):
            if st.startswith(("switch", "case", "otherwise", "if ")):
                skipped.append(st[:80])
            continue
        m = _LHS.match(st)
        if not m:
            skipped.append(st[:80])
            continue
        name, idx, rhs = m.group(1), m.group(3), st[m.end():]
        if re.match(r"\s*setup_\w+\(MWS\)", rhs) or "'" in rhs and re.search(r"'[A-Za-z]", rhs):
            skipped.append(st[:80])          # calls to other setup files, or strings
            continue
        try:
            p = _Parser(_tokenize(rhs), env)
            val = p.expr()
            if p.peek() is not None:
                raise ValueError("trailing tokens: %r" % p.peek().text)
        except Exception as exc:           # recorded, never silently ignored
            skipped.append("%s  <-- %s" % (st[:80].replace("\n", " "), exc))
            continue
        if idx is None:
            env[name] = val
            continue
        parts = [s.strip() for s in idx.split(",")]
        arr = env.get(name)
        v = np.atleast_2d(np.asarray(val, dtype=float))
        if len(parts) == 2 and parts[1] == ":":            # name(i,:) = row
            i = int(parts[0]) - 1
            if arr is None or not isinstance(arr, np.ndarray):
                arr = np.zeros((0, v.size))
            if arr.shape[0] <= i:
                arr = np.vstack([arr, np.zeros((i + 1 - arr.shape[0], v.size))])
            arr[i, :] = v.reshape(-1)
        elif len(parts) == 2 and parts[0] == ":":          # name(:,j) = column
            j = int(parts[1]) - 1
            if arr is None or not isinstance(arr, np.ndarray):
                arr = np.zeros((v.size, 0))
            if arr.shape[1] <= j:
                arr = np.hstack([arr, np.zeros((v.size, j + 1 - arr.shape[1]))])
            arr[:, j] = v.reshape(-1)
        elif len(parts) == 1:                              # name(k) = scalar
            k = int(parts[0]) - 1
            if arr is None or not isinstance(arr, np.ndarray):
                arr = np.zeros((1, 0))
            if arr.shape[1] <= k:
                arr = np.hstack([arr, np.zeros((1, k + 1 - arr.shape[1]))])
            arr[0, k] = float(v.reshape(-1)[0])
        else:
            skipped.append(st[:80])
            continue
        env[name] = arr


# --------------------------------------------------------------------------
# 2. T-MATS gas-table reader (C arrays)
# --------------------------------------------------------------------------
_C_ARRAY = re.compile(r"(?:double|int)\s+(\w+)\[(\d+)\]\s*=\s*\{([^}]*)\}", re.S)


def read_gas_tables(mex_dir: Path) -> dict:
    """Collect the polynomial coefficient arrays of t2hc/h2tc/pt2sc/sp2tc and check they agree between files."""
    tables: dict[str, list[float]] = {}
    for fname in ("t2hc_TMATS.c", "h2tc_TMATS.c", "pt2sc_TMATS.c", "sp2tc_TMATS.c"):
        src = (mex_dir / fname).read_text(encoding="latin-1")
        for name, n, body in _C_ARRAY.findall(src):
            vals = [float(x) for x in re.split(r"[,\s]+", body.strip()) if x]
            if len(vals) != int(n):
                raise ValueError("%s: array %s has %d values, expected %s" % (fname, name, len(vals), n))
            if name in tables and tables[name] != vals:
                raise ValueError("gas table %s differs between T-MATS files" % name)
            tables[name] = vals
    return tables


def read_ambient_tables(eng_mdl: Path) -> dict:
    """Standard-atmosphere tables of the Ambient block, read from its mask in AGTF30_eng.mdl.

    The .mdl is a text file (CR line endings); long strings continue on the next
    line as a second quoted string.
    """
    text = eng_mdl.read_bytes().decode("latin-1").replace("\r\n", "\n").replace("\r", "\n")
    out = {}
    for key, name in (("X_A_AltVec_M", "AltVec"), ("T_A_TsVec_M", "TsVec"), ("T_A_PsVec_M", "PsVec")):
        m = re.search(key + r'\s+((?:"[^"]*"\s*)+)', text)
        if not m:
            raise ValueError("Ambient table %s not found in %s" % (key, eng_mdl))
        expr = "".join(re.findall(r'"([^"]*)"', m.group(1)))
        val = _Parser(_tokenize(expr), {}).expr()
        out[name] = np.asarray(val, dtype=float).reshape(-1).tolist()
    return out


# --------------------------------------------------------------------------
# 3. Public loader with a JSON cache
# --------------------------------------------------------------------------
def parse_vendor(vendor_dir: Path = VENDOR) -> dict:
    """Parse everything from the vendored NASA sources. Returns {'MWS': {...}, 'gas': {...}, 'skipped': [...]}."""
    env: dict = {}
    skipped: list[str] = []
    setup = vendor_dir / "agtf30" / "AGTF30" / "SimSetup"
    for fname in SETUP_FILES:
        text = (setup / fname).read_bytes().decode("latin-1")
        run_matlab(text, env, skipped)
    mws = {}
    for key, val in env.items():
        if not key.startswith("MWS."):
            continue
        if isinstance(val, np.ndarray):
            val = val.tolist()
        mws[key[4:]] = val
    gas = read_gas_tables(vendor_dir / "tmats" / "TMATS_Library" / "MEX")
    ambient = read_ambient_tables(vendor_dir / "agtf30" / "AGTF30" / "AGTF30_eng.mdl")
    return {"MWS": mws, "gas": gas, "ambient": ambient, "skipped": skipped}


def load_params(vendor_dir: Path | str | None = None, use_cache: bool = True) -> dict:
    """Return the parsed NASA parameters; cached as JSON in eigenmap/data (re-created from vendor if missing)."""
    vendor = Path(vendor_dir) if vendor_dir else VENDOR
    cache = DATA_DIR / "agtf30_params.json"
    if use_cache and vendor == VENDOR and cache.exists():
        return json.loads(cache.read_text())
    params = parse_vendor(vendor)
    if vendor == VENDOR:
        DATA_DIR.mkdir(parents=True, exist_ok=True)
        cache.write_text(json.dumps(params, indent=0, sort_keys=True))
    return params


if __name__ == "__main__":
    prm = load_params(use_cache=False)
    print("parsed %d MWS entries, %d gas tables, %d skipped statements" % (len(prm["MWS"]), len(prm["gas"]), len(prm["skipped"])))
    for s in prm["skipped"]:
        print("  skipped:", s.replace("\n", " "))
