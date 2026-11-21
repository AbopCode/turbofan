"""Mini evaluator for the small MATLAB subset used by the AGTF30 setup_*.m files.

Only what those files (and the Simulink mask strings) need: numbers, matrix
literals, ranges a:s:b, + - * / ^, transpose, struct fields, simple indexing with
':' and 'end', and the functions ones/zeros/max/min/length/numel.
Everything is held as 2-D (or 3-D) numpy arrays, like MATLAB does.
"""
from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Any, Callable

import numpy as np


class MatlabError(Exception):
    """Raised when a statement is outside the supported subset."""


def read_text(path: str | Path) -> str:
    """Read a text file that may use CR, LF or CRLF line endings; return LF text."""
    with open(path, "r", newline="", encoding="latin-1") as fh:
        raw = fh.read()
    return re.sub(r"\r\n|\r", "\n", raw)


# --------------------------------------------------------------------------- tokens
_NUM = re.compile(r"(\d+\.?\d*([eE][+-]?\d+)?|\.\d+([eE][+-]?\d+)?)")
_ID = re.compile(r"[A-Za-z_]\w*(\.[A-Za-z_]\w*)*")
_OPERAND_END = ("num", "id", ")", "]", "'")


def tokenize(src: str) -> list[tuple[str, Any]]:
    """Split an expression into tokens; whitespace separates elements inside [ ]."""
    toks: list[tuple[str, Any]] = []
    stack: list[str] = []
    i, n = 0, len(src)
    while i < n:
        c = src[i]
        if c in " \t":
            j = i
            while j < n and src[j] in " \t":
                j += 1
            in_br = bool(stack) and stack[-1] == "["
            if in_br and toks and toks[-1][0] in _OPERAND_END and j < n:
                d = src[j]
                if d in "+-":
                    if j + 1 < n and src[j + 1] not in " \t":
                        toks.append((",", None))
                elif d.isalnum() or d in "._([":
                    toks.append((",", None))
            i = j
            continue
        if c == "\n":
            if stack and stack[-1] == "[":
                toks.append((";", None))
            i += 1
            continue
        m = _NUM.match(src, i)
        if m and (c.isdigit() or c == "."):
            toks.append(("num", float(m.group(0))))
            i = m.end()
            continue
        m = _ID.match(src, i)
        if m:
            toks.append(("id", m.group(0)))
            i = m.end()
            continue
        if c in "([":
            stack.append(c)
            toks.append((c, None))
        elif c in ")]":
            if not stack:
                raise MatlabError("unbalanced bracket in: " + src[:80])
            stack.pop()
            toks.append((c, None))
        elif c == "'":
            if toks and toks[-1][0] in _OPERAND_END:
                toks.append(("'", None))
            else:
                raise MatlabError("strings are not supported: " + src[:80])
        elif c == "." and i + 1 < n and src[i + 1] in "*/^":
            toks.append((src[i + 1], None))
            i += 1
        elif c in "+-*/^:,;":
            toks.append((c, None))
        else:
            raise MatlabError(f"unexpected character {c!r} in: {src[:80]}")
        i += 1
    return toks


# --------------------------------------------------------------------------- values
def _arr(v: Any) -> np.ndarray:
    a = np.asarray(v, dtype=float)
    if a.ndim == 0:
        a = a.reshape(1, 1)
    elif a.ndim == 1:
        a = a.reshape(1, -1)
    return a


def _colon(a: float, s: float, b: float) -> np.ndarray:
    """MATLAB colon a:s:b (values built from both ends, as MATLAB does)."""
    if s == 0 or (b - a) / s < -1e-12:
        return np.zeros((1, 0))
    n = int(math.floor((b - a) / s + 1e-10)) + 1
    last = a + (n - 1) * s
    out = [0.0] * n
    for k in range(n):
        out[k] = a + k * s if k < (n + 1) // 2 else last - (n - 1 - k) * s
    return _arr(out)


class Evaluator:
    """Evaluates expressions / runs statements against a workspace dict."""

    def __init__(self, ws: dict[str, Any], call_hook: Callable[[str, list], Any] | None = None):
        self.ws = ws
        self.call_hook = call_hook
        self.toks: list[tuple[str, Any]] = []
        self.p = 0
        self._end: list[int] = []

    # -- workspace access
    def lookup(self, name: str) -> Any:
        parts = name.split(".")
        if parts[0] not in self.ws:
            raise MatlabError("unknown variable " + name)
        v = self.ws[parts[0]]
        for k in parts[1:]:
            if not isinstance(v, dict) or k not in v:
                raise MatlabError("unknown field " + name)
            v = v[k]
        return v

    def has(self, name: str) -> bool:
        try:
            self.lookup(name)
            return True
        except MatlabError:
            return False

    def store(self, name: str, value: Any) -> None:
        parts = name.split(".")
        d = self.ws
        for k in parts[:-1]:
            if k not in d or not isinstance(d[k], dict):
                d[k] = {}
            d = d[k]
        d[parts[-1]] = value

    # -- parser helpers
    def peek(self) -> str | None:
        return self.toks[self.p][0] if self.p < len(self.toks) else None

    def take(self, kind: str | None = None) -> tuple[str, Any]:
        if self.p >= len(self.toks):
            raise MatlabError("unexpected end of expression")
        t = self.toks[self.p]
        if kind is not None and t[0] != kind:
            raise MatlabError(f"expected {kind!r}, got {t[0]!r}")
        self.p += 1
        return t

    def eval(self, src: str) -> np.ndarray:
        """Evaluate one expression string and return a numpy array (>= 2-D)."""
        self.toks = tokenize(src)
        self.p = 0
        v = self.expr()
        if self.p != len(self.toks):
            raise MatlabError("trailing tokens in: " + src[:80])
        return v

    def expr(self) -> np.ndarray:
        a = self.additive()
        if self.peek() == ":":
            self.take()
            b = self.additive()
            if self.peek() == ":":
                self.take()
                c = self.additive()
                return _colon(float(a[0, 0]), float(b[0, 0]), float(c[0, 0]))
            return _colon(float(a[0, 0]), 1.0, float(b[0, 0]))
        return a

    def additive(self) -> np.ndarray:
        a = self.term()
        while self.peek() in ("+", "-"):
            op = self.take()[0]
            b = self.term()
            a = a + b if op == "+" else a - b
        return a

    def term(self) -> np.ndarray:
        a = self.unary()
        while self.peek() in ("*", "/"):
            op = self.take()[0]
            b = self.unary()
            if a.size != 1 and b.size != 1 and a.shape != b.shape:
                raise MatlabError("matrix product not supported")
            a = a * b if op == "*" else a / b
        return a

    def unary(self) -> np.ndarray:
        if self.peek() == "-":
            self.take()
            return -self.unary()
        if self.peek() == "+":
            self.take()
            return self.unary()
        return self.power()

    def power(self) -> np.ndarray:
        a = self.postfix()
        if self.peek() == "^":
            self.take()
            sign = 1.0
            while self.peek() in ("-", "+"):
                if self.take()[0] == "-":
                    sign = -sign
            b = self.postfix()
            return _arr(np.power(a, sign * b))
        return a

    def postfix(self) -> np.ndarray:
        a = self.primary()
        while self.peek() == "'":
            self.take()
            a = a.T.copy()
        return a

    def primary(self) -> np.ndarray:
        k, v = self.take()
        if k == "num":
            return _arr(v)
        if k == "(":
            a = self.expr()
            self.take(")")
            return a
        if k == "[":
            return self.matrix()
        if k == "id":
            if v == "end":
                if not self._end:
                    raise MatlabError("'end' outside an index")
                return _arr(float(self._end[-1]))
            if v == "pi":
                return _arr(math.pi)
            if self.peek() == "(":
                return self.call_or_index(v)
            return _arr(self.lookup(v))
        raise MatlabError(f"unexpected token {k!r}")

    def matrix(self) -> np.ndarray:
        rows: list[list[np.ndarray]] = [[]]
        while True:
            k = self.peek()
            if k == "]":
                self.take()
                break
            if k == ";":
                self.take()
                if rows[-1]:
                    rows.append([])
                continue
            if k == ",":
                self.take()
                continue
            rows[-1].append(self.expr())
        rows = [r for r in rows if r]
        if not rows:
            return np.zeros((0, 0))
        return np.vstack([np.hstack(r) for r in rows])

    def call_or_index(self, name: str) -> np.ndarray:
        self.take("(")
        is_var = self.has(name)
        base = _arr(self.lookup(name)) if is_var else None
        args: list[Any] = []
        idx = 0
        while self.peek() != ")":
            if self.peek() == ",":
                self.take()
                idx += 1
                continue
            nxt = self.toks[self.p + 1][0] if self.p + 1 < len(self.toks) else None
            if self.peek() == ":" and nxt in (",", ")"):
                self.take()
                args.append(":")
                continue
            if is_var:
                self._end.append(self._end_value(base, idx))
                args.append(self.expr())
                self._end.pop()
            else:
                args.append(self.expr())
        self.take(")")
        if is_var:
            return self._index(base, args)
        return self._call(name, args)

    def _end_value(self, base: np.ndarray, idx: int) -> int:
        # one index -> numel; two indices -> size along that dimension (heuristic: look ahead)
        depth, n_args, q = 0, 1, self.p
        while q < len(self.toks):
            k = self.toks[q][0]
            if k in "([":
                depth += 1
            elif k in ")]":
                if depth == 0:
                    break
                depth -= 1
            elif k == "," and depth == 0:
                n_args += 1
            q += 1
        total = idx + n_args
        return int(base.size) if total == 1 else int(base.shape[idx])

    @staticmethod
    def _index(base: np.ndarray, args: list[Any]) -> np.ndarray:
        def sel(a: Any, n: int) -> np.ndarray:
            if isinstance(a, str):
                return np.arange(n)
            return np.asarray(a, dtype=float).ravel().round().astype(int) - 1

        if len(args) == 1:
            flat = base.ravel(order="F")
            out = flat[sel(args[0], flat.size)]
            return out.reshape(-1, 1) if (base.shape[1] == 1 and base.shape[0] > 1) else out.reshape(1, -1)
        if len(args) == 2:
            return base[np.ix_(sel(args[0], base.shape[0]), sel(args[1], base.shape[1]))]
        raise MatlabError("indexing with more than 2 subscripts is not supported")

    def _call(self, name: str, args: list[Any]) -> np.ndarray:
        if name in ("ones", "zeros"):
            r = int(args[0][0, 0])
            c = int(args[1][0, 0]) if len(args) > 1 else r
            return np.ones((r, c)) if name == "ones" else np.zeros((r, c))
        if name in ("max", "min"):
            f = np.max if name == "max" else np.min
            if len(args) == 2:
                return _arr(np.maximum(args[0], args[1]) if name == "max" else np.minimum(args[0], args[1]))
            a = args[0]
            if a.shape[0] == 1 or a.shape[1] == 1:
                return _arr(f(a))
            return _arr(f(a, axis=0))
        if name == "length":
            return _arr(float(max(args[0].shape) if args[0].size else 0))
        if name == "numel":
            return _arr(float(args[0].size))
        if self.call_hook is not None:
            return self.call_hook(name, args)
        raise MatlabError("unknown function " + name)

    # -- statements
    def assign(self, lhs: str, index: str | None, value: np.ndarray) -> None:
        """Assign to a (possibly indexed) variable, growing it like MATLAB does."""
        if index is None:
            self.store(lhs, value)
            return
        subs = [s.strip() for s in index.split(",")]
        cur = _arr(self.lookup(lhs)) if self.has(lhs) else np.zeros((0, 0))

        def sel(s: str, n_have: int, n_val: int) -> np.ndarray:
            if s == ":":
                return np.arange(n_have if n_have else n_val)
            return np.array([int(round(float(s))) - 1])

        if len(subs) == 2:
            r = sel(subs[0], cur.shape[0], value.shape[0])
            c = sel(subs[1], cur.shape[1], value.shape[1])
            nr, nc = max(cur.shape[0], r.max() + 1), max(cur.shape[1], c.max() + 1)
            if (nr, nc) != cur.shape:
                new = np.zeros((nr, nc))
                new[: cur.shape[0], : cur.shape[1]] = cur
                cur = new
            cur[np.ix_(r, c)] = value.reshape(len(r), len(c))
        elif len(subs) == 3 and subs[0] == ":" and subs[1] == ":":
            k = int(subs[2]) - 1
            if cur.ndim != 3:
                cur = np.zeros(value.shape + (0,))
            if cur.shape[2] <= k:
                new = np.zeros(value.shape + (k + 1,))
                new[:, :, : cur.shape[2]] = cur
                cur = new
            cur[:, :, k] = value
        else:
            raise MatlabError("unsupported indexed assignment " + lhs + "(" + index + ")")
        self.store(lhs, cur)


_ASSIGN = re.compile(r"^([A-Za-z_]\w*(?:\.[A-Za-z_]\w*)*)\s*(?:\(([^=]*)\))?\s*=(?!=)\s*(.*)$", re.S)


def split_statements(text: str) -> list[str]:
    """Strip comments / continuations and split a script into statements."""
    lines = []
    for ln in text.split("\n"):
        k = ln.find("%")
        if k >= 0:
            ln = ln[:k]
        lines.append(ln.rstrip())
    src = "\n".join(lines)
    src = re.sub(r"\.\.\.[^\n]*\n", " ", src)       # line continuation
    out, buf, depth = [], [], 0
    for ch in src:
        if ch in "([":
            depth += 1
        elif ch in ")]":
            depth -= 1
        if depth == 0 and ch in ";\n":
            s = "".join(buf).strip()
            if s:
                out.append(s)
            buf = []
        else:
            buf.append(ch)
    s = "".join(buf).strip()
    if s:
        out.append(s)
    return out


def run_script(path: str | Path, ws: dict[str, Any], strict: bool = True,
               stop_at: str | None = None) -> list[str]:
    """Run a setup_*.m script into workspace `ws` (MWS lives in ws['MWS']).

    Local variables of the script are kept in a private scope; `MWS = setup_X(MWS)`
    runs setup_X.m from the same folder. Returns the list of skipped statements.
    """
    path = Path(path)
    text = read_text(path)
    if stop_at is not None and stop_at in text:
        text = text[: text.index(stop_at)]
    local: dict[str, Any] = {"MWS": ws.setdefault("MWS", {})}
    skipped: list[str] = []
    ev = Evaluator(local)
    for st in split_statements(text):
        if st.startswith("function ") or st == "end":
            continue
        m = re.match(r"^MWS\s*=\s*(setup_\w+)\s*\(\s*MWS\s*\)$", st)
        if m:
            sub = path.with_name(m.group(1) + ".m")
            if sub.name in ("setup_Sensor.m", "setup_Actuators.m"):
                continue                                  # not part of the engine physics
            skipped += run_script(sub, ws, strict=strict)
            continue
        m = _ASSIGN.match(st)
        try:
            if not m:
                raise MatlabError("not an assignment: " + st[:60])
            ev.assign(m.group(1), m.group(2), ev.eval(m.group(3)))
        except MatlabError as exc:
            if strict:
                raise MatlabError(f"{path.name}: {exc} :: {st[:80]}") from exc
            skipped.append(st[:80])
    return skipped


def to_python(v: Any) -> Any:
    """Convert workspace values: 1x1 -> float, vectors -> 1-D array, dict recursively."""
    if isinstance(v, dict):
        return {k: to_python(x) for k, x in v.items()}
    a = np.asarray(v, dtype=float)
    if a.size == 1:
        return float(a.ravel()[0])
    if a.ndim == 2 and 1 in a.shape:
        return a.ravel().copy()
    return a.copy()
