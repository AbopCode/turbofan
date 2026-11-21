"""Export outputs/engine_viewer.html: a standalone, offline replay of the animated engine (no Python, no internet).

It records a throttle burst (sea level, idle -> max at t = 1 s, NASA baseline controller) with the lab's live simulation
and embeds the frames plus the lab's own canvas animation (twin/assets/engine_canvas.js) in one HTML file.
Then, if Chrome/Chromium is available, it takes screenshots of the lab (Simulator, Live engine) and of the viewer.
EXTENSION - outside the report objectives.

Run: python -m twin.export_viewer
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

from eigenmap import config

from . import sim

HERE = Path(__file__).resolve().parent

PAGE = """<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><title>AGTF30 engine viewer (offline replay)</title>
<style>
body{font-family:system-ui,-apple-system,sans-serif;margin:0;background:#f2f2ef;color:#161616}
header{background:#14213d;color:#fff;padding:10px 18px} header b{font-size:18px} header span{color:#c9d3e6;margin-left:12px;font-size:13px}
main{padding:14px 18px} #engine-canvas{width:100%;height:520px;display:block;border:1px solid #dcdcd6;border-radius:8px;background:#f6f7f8}
.bar{display:flex;gap:10px;align-items:center;margin:10px 0} button{font-weight:600;padding:7px 14px;border-radius:6px;border:1px solid #c9c8c2;background:#fff;cursor:pointer}
button.go{background:#1baf7a;color:#fff;border-color:#138a5f} input[type=range]{flex:1} .muted{color:#666;font-size:12.5px}
</style></head><body>
<header><b>AGTF30 engine viewer</b><span>offline replay of a sea-level throttle burst (idle &rarr; max at t = 1 s, NASA baseline controller) &middot; schematic, not an engineering drawing</span></header>
<main>
<canvas id="engine-canvas"></canvas>
<div class="bar"><button class="go" id="play">&#9654; Play</button><button id="pause">&#10073;&#10073; Pause</button>
<input type="range" id="pos" min="0" max="__N__" value="0"><span id="time">t = 0.0 s</span>
<label class="muted"><input type="checkbox" id="labels" checked> labels</label></div>
<p class="muted">Recorded with the Turbofan Virtual Lab (twin/). For live experiments run <code>./run.sh app</code>. Rotation is shown about 30x slower than real.</p>
</main>
<script>__JS__</script>
<script>
const FR = __FRAMES__; let i = 0, playing = false;
const pos = document.getElementById('pos'), tl = document.getElementById('time');
function show(k){ const f = Object.assign({}, FR[k], {running: true}); LabEngine.update(f); pos.value = k; tl.textContent = 't = ' + f.t.toFixed(1) + ' s'; }
document.getElementById('play').onclick = () => { playing = true; if (i >= FR.length - 1) i = 0; };
document.getElementById('pause').onclick = () => { playing = false; };
pos.oninput = () => { i = +pos.value; show(i); };
document.getElementById('labels').onchange = (e) => LabEngine.setLabels(e.target.checked);
setInterval(() => { if (playing) { i = Math.min(i + 1, FR.length - 1); show(i); if (i >= FR.length - 1) playing = false; } }, 100);
show(0); playing = true;
</script></body></html>"""


def record(t_end: float = 14.0) -> list:
    """Frames every 0.1 s of a sea-level burst, recorded with the lab's live simulation."""
    L = sim.LiveSim()
    L.update({"mode": "auto", "power": 0.0, "alt": 0.0, "mach": 0.0, "dT": 0.0})
    L.reset()
    L.running = True
    frames = []
    n = int(round(t_end / sim.DT))
    for k in range(n + 1):
        if abs(L.t - 1.0) < sim.DT / 2:
            L.update({"power": 100.0})
        if k % 5 == 0:
            snap = L.snapshot()
            frames.append({key: (round(v, 5) if isinstance(v, float) else v) for key, v in snap.items()})
        if not L._one_step():
            break
    return frames


def main() -> int:
    frames = record()
    js = (HERE / "assets" / "engine_canvas.js").read_text()
    html = PAGE.replace("__JS__", js).replace("__FRAMES__", json.dumps(frames, separators=(",", ":"))).replace("__N__", str(len(frames) - 1))
    config.OUT.mkdir(parents=True, exist_ok=True)
    path = config.OUT / "engine_viewer.html"
    path.write_text(html)
    print("PASS export: wrote engine_viewer.html (%d frames, %.0f kB, works offline, no WebGL)" % (len(frames), path.stat().st_size / 1e3))
    from .screenshot import app_screenshot, find_chrome, shoot
    figs = config.FIGS
    figs.mkdir(parents=True, exist_ok=True)
    shots = [("?tab=live&scenario=burst", figs / "lab_live_engine.png", "1600,1500", 30000),
             ("?preset=burst", figs / "lab_simulator_burst.png", "1600,1900", 60000),
             ("?tab=eigmap&run=1", figs / "lab_eigenvalue_map.png", "1600,1500", 60000)]
    if not find_chrome():
        print("SKIP export: no Chrome/Chromium found - lab screenshots not taken")
        return 0
    served, ok = app_screenshot(shots, config.ROOT)
    print("%s export: the virtual lab launched headlessly and served its page (HTTP 200)" % ("PASS" if served else "FAIL"))
    for (q, out, _s, _b), good in zip(shots, ok):
        print("%s export: screenshot %s (%s)" % ("PASS" if good else "FAIL", out.name, q))
    ok_v = shoot(path.resolve().as_uri(), figs / "lab_offline_viewer.png", size="1500,800", budget_ms=8000)
    print("%s export: screenshot lab_offline_viewer.png" % ("PASS" if ok_v else "FAIL"))
    return 0 if (served and all(ok) and ok_v) else 1


if __name__ == "__main__":
    sys.exit(main())
