#!/usr/bin/env bash
# ============================================================================
#  The Eigenvalue Map of a Turbofan Engine - one-command runner (macOS / Linux)
#
#    ./run.sh                 quick verify, then start the virtual lab and open the browser
#    ./run.sh verify          quick check: tests + PASS/FAIL table            (< 60 s)
#    ./run.sh full [--fast]   run the whole study, then the table  ([--jobs N] workers)
#    ./run.sh app             start the virtual lab only (http://127.0.0.1:8050)
#    ./run.sh export          write outputs/engine_viewer.html (offline engine replay) + lab screenshots
#    ./run.sh fetch-data      optional: download NASA N-CMAPSS data into data/ ([--yes])
#    ./run.sh clean           delete generated files (never source, vendor/ or .venv)
#
#  Environment variables (all optional):
#    EIGENMAP_PYTHON=/path/to/python   use this interpreter, do not create .venv
#    PORT=8050                         port of the virtual lab
#    EIGENMAP_NO_BROWSER=1             do not open a browser window
#    EIGENMAP_OUT=/some/folder         write outputs there instead of ./outputs
#
#  Kept compatible with the old bash 3.2 that ships with macOS.
# ============================================================================
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"          # works when called from any directory
ROOT="$(pwd)"
export MPLBACKEND=Agg                         # matplotlib without a display

VENV="$ROOT/.venv"
STAMP="$VENV/.requirements.sha256"            # hash of requirements.txt at the last install
PORT="${PORT:-8050}"
PY=""                                         # interpreter used for the project (set below)

# ---- coloured status lines (colours only on a real terminal) ----------------
if [ -t 1 ] && [ -z "${NO_COLOR:-}" ]; then
  C_RED=$'\033[31m'; C_GRN=$'\033[32m'; C_YEL=$'\033[33m'; C_CYN=$'\033[36m'; C_OFF=$'\033[0m'
else
  C_RED=""; C_GRN=""; C_YEL=""; C_CYN=""; C_OFF=""
fi
info() { printf '%s[run]%s %s\n'  "$C_CYN" "$C_OFF" "$*"; }
ok()   { printf '%s[ ok]%s %s\n'  "$C_GRN" "$C_OFF" "$*"; }
warn() { printf '%s[warn]%s %s\n' "$C_YEL" "$C_OFF" "$*" >&2; }
fail() { printf '%s[fail]%s %s\n' "$C_RED" "$C_OFF" "$*" >&2; }

usage() {
  # print the comment block at the top of this file
  sed -n '2,20p' "$0" | sed 's/^#//'
}

# ---- small helpers ------------------------------------------------------------
file_hash() {                                  # sha256 of a file (Linux, macOS, or via python)
  if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1" | cut -d' ' -f1
  elif command -v shasum >/dev/null 2>&1; then shasum -a 256 "$1" | cut -d' ' -f1
  else "$PY" -c 'import hashlib,sys; print(hashlib.sha256(open(sys.argv[1],"rb").read()).hexdigest())' "$1"
  fi
}

find_python() {                                # print the first python of the wanted version
  local version="$1" candidate
  for candidate in "python$version" "/opt/homebrew/bin/python$version" "/usr/local/bin/python$version"; do
    if command -v "$candidate" >/dev/null 2>&1; then command -v "$candidate"; return 0; fi
  done
  return 1
}

# ---- environment bootstrap: create .venv once, install requirements when they change ----
create_venv() {
  local base=""
  if command -v uv >/dev/null 2>&1; then
    info "creating .venv with uv (Python 3.13 preferred, 3.12 as fallback)"
    if uv venv --python 3.13 "$VENV" || uv venv --python 3.12 "$VENV"; then return 0; fi
    warn "uv could not create the environment; trying the Python installations on this machine"
    rm -rf "$VENV"
  fi
  if base="$(find_python 3.13)" || base="$(find_python 3.12)"; then
    info "creating .venv with $base"
  elif command -v python3 >/dev/null 2>&1; then
    base="$(command -v python3)"
    warn "Python 3.13 / 3.12 not found. Falling back to the system python3:"
    warn "  $base ($("$base" --version 2>&1)). The study is tested on 3.13 and 3.12 only."
    if ! "$base" -c 'import sys; sys.exit(0 if sys.version_info >= (3, 10) else 1)'; then
      fail "this python3 is older than 3.10 and cannot run the project."
      fail "Install Python 3.13 (macOS: 'brew install python@3.13', or https://www.python.org/downloads/)."
      exit 1
    fi
  else
    fail "no Python found. Install Python 3.13 (macOS: 'brew install python@3.13')."
    exit 1
  fi
  "$base" -m venv "$VENV"
}

install_requirements() {
  local want have=""
  want="$(file_hash "$ROOT/requirements.txt")"
  if [ -f "$STAMP" ]; then have="$(cat "$STAMP")"; fi
  if [ "$want" = "$have" ]; then return 0; fi            # nothing changed: skip the install
  info "installing Python packages from requirements.txt (first run, or the file changed)"
  if command -v uv >/dev/null 2>&1; then
    uv pip install --python "$PY" -r "$ROOT/requirements.txt"
  else
    # an environment made by uv has no pip of its own: add it if needed
    "$PY" -m pip --version >/dev/null 2>&1 || "$PY" -m ensurepip --upgrade
    "$PY" -m pip install --disable-pip-version-check -r "$ROOT/requirements.txt"
  fi
  printf '%s\n' "$want" > "$STAMP"
  ok "packages installed"
}

bootstrap() {
  if [ -n "${EIGENMAP_PYTHON:-}" ]; then                  # explicit interpreter: no .venv
    PY="$EIGENMAP_PYTHON"
    info "using EIGENMAP_PYTHON=$PY (no virtual environment is created)"
    return 0
  fi
  PY="$VENV/bin/python"
  if [ ! -x "$PY" ]; then create_venv; fi
  if [ ! -x "$PY" ]; then fail "could not create $VENV"; exit 1; fi
  install_requirements
}

# ---- commands -------------------------------------------------------------------
wait_then_open_browser() {                     # runs in the background next to the app
  local url="$1" i
  for i in $(seq 1 60); do                     # wait (max 30 s) until the port answers
    if "$PY" -c 'import socket,sys; socket.create_connection(("127.0.0.1", int(sys.argv[1])), 0.5).close()' "$PORT" 2>/dev/null; then
      break
    fi
    sleep 0.5
  done
  case "$(uname -s)" in
    Darwin) open "$url" >/dev/null 2>&1 || true ;;
    *)      if command -v xdg-open >/dev/null 2>&1; then xdg-open "$url" >/dev/null 2>&1 || true; fi ;;
  esac
}

run_app() {
  local url="http://127.0.0.1:$PORT/"
  if [ ! -f "$ROOT/twin/app.py" ]; then
    warn "the virtual lab (twin/app.py) is missing - nothing to launch."
    return 0
  fi
  info "starting the virtual lab at $url   (stop it with Ctrl-C)"
  if [ -z "${EIGENMAP_NO_BROWSER:-}" ]; then wait_then_open_browser "$url" & fi
  exec "$PY" -m twin.app --port "$PORT"        # the app replaces this shell, so Ctrl-C reaches it
}

fetch_data() {
  # OPTIONAL. The study does not need this data; every failure here is non-fatal (exit 0).
  local url='https://phm-datasets.s3.amazonaws.com/NASA/17.+Turbofan+Engine+Degradation+Simulation+Data+Set+2.zip'
  local page='https://www.nasa.gov/intelligent-systems-division/discovery-and-systems-health/pcoe/pcoe-data-set-repository/'
  local target="$ROOT/data/N-CMAPSS_Turbofan_Engine_Degradation_Simulation_Data_Set_2.zip"
  local yes="no" answer=""
  if [ "${1:-}" = "--yes" ] || [ "${1:-}" = "-y" ]; then yes="yes"; fi

  info "NASA N-CMAPSS: 'Turbofan Engine Degradation Simulation Data Set-2' (NASA Prognostics Data Repository)"
  warn "the download address below could NOT be verified when this script was written:"
  echo "         $url"
  echo "       If it fails, download by hand from the repository page (data set no. 17):"
  echo "         $page"
  echo "       and put the zip file into: $ROOT/data/"
  warn "the file is LARGE (several GB). The study and the lab do NOT need it."

  if [ "$yes" != "yes" ]; then
    if [ ! -t 0 ]; then
      warn "not an interactive terminal and no --yes given: download skipped."
      return 0
    fi
    printf 'Download now? [y/N] '
    read -r answer || answer=""
    case "$answer" in y|Y|yes|YES) ;; *) info "download skipped."; return 0 ;; esac
  fi

  mkdir -p "$ROOT/data"
  if command -v curl >/dev/null 2>&1; then
    # -f: fail on HTTP errors, -L: follow redirects, -C -: resume a partial download
    if curl -fL -C - --retry 2 -o "$target.part" "$url"; then mv "$target.part" "$target"; fi
  elif command -v wget >/dev/null 2>&1; then
    if wget -c -O "$target.part" "$url"; then mv "$target.part" "$target"; fi
  else
    warn "neither curl nor wget is installed."
  fi
  if [ -f "$target" ]; then
    ok "saved $target"
    info "unpack it with:  unzip '$target' -d '$ROOT/data'   (reading the .h5 files needs: pip install h5py)"
  else
    warn "download failed (this is not an error for the study). Manual download page:"
    echo "         $page"
  fi
  return 0
}

clean() {
  # Needs no packages, so do not build .venv just to clean.
  local py=""
  if [ -n "${EIGENMAP_PYTHON:-}" ]; then py="$EIGENMAP_PYTHON"
  elif [ -x "$VENV/bin/python" ]; then py="$VENV/bin/python"
  elif command -v python3 >/dev/null 2>&1; then py="$(command -v python3)"
  else fail "no Python found to run the clean command."; exit 1
  fi
  "$py" -m eigenmap.cli clean
}

command="${1:-default}"
if [ $# -gt 0 ]; then shift; fi

case "$command" in
  default)
    bootstrap
    if ! "$PY" -m eigenmap.cli verify; then
      fail "verification failed (see the table above) - the lab was not started."
      fail "Start it anyway with:  ./run.sh app"
      exit 1
    fi
    ok "verification passed"
    run_app
    ;;
  verify)      bootstrap; "$PY" -m eigenmap.cli verify "$@" ;;
  full)        bootstrap; "$PY" -m eigenmap.cli full "$@" ;;
  app)         bootstrap; run_app ;;
  export)      bootstrap; "$PY" -m eigenmap.cli export "$@" ;;
  fetch-data)  fetch_data "$@" ;;
  clean)       clean ;;
  help|-h|--help) usage ;;
  *)           fail "unknown command: $command"; usage; exit 2 ;;
esac
