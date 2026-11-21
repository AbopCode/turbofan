@echo off
rem ============================================================================
rem  The Eigenvalue Map of a Turbofan Engine - one-command runner (Windows)
rem
rem    run.bat                 quick verify, then start the 3D twin and open the browser
rem    run.bat verify          quick check: tests + PASS/FAIL table            (< 60 s)
rem    run.bat full [--fast]   run the whole study, then the table  ([--jobs N] workers)
rem    run.bat app             start the twin only (http://127.0.0.1:8050)
rem    run.bat export          write the stand-alone viewer outputs\engine_viewer.html
rem    run.bat fetch-data      optional: download NASA N-CMAPSS data into data\ ([--yes])
rem    run.bat clean           delete generated files (never source, vendor\ or .venv)
rem
rem  Optional environment variables: PORT (default 8050), EIGENMAP_NO_BROWSER=1,
rem  EIGENMAP_OUT=folder for the outputs.
rem ============================================================================
setlocal EnableExtensions
cd /d "%~dp0"
set "MPLBACKEND=Agg"
set "VENV=%CD%\.venv"
set "PY=%VENV%\Scripts\python.exe"
set "STAMP=%VENV%\.requirements.sha256"
if not defined PORT set "PORT=8050"

set "CMD=%~1"
if "%CMD%"=="" set "CMD=default"

rem Everything after the command name is passed on (at most 8 extra words).
set "ARGS=%2 %3 %4 %5 %6 %7 %8 %9"

if /i "%CMD%"=="default"    goto :cmd_default
if /i "%CMD%"=="verify"     goto :cmd_verify
if /i "%CMD%"=="full"       goto :cmd_full
if /i "%CMD%"=="app"        goto :cmd_app
if /i "%CMD%"=="export"     goto :cmd_export
if /i "%CMD%"=="fetch-data" goto :cmd_fetch
if /i "%CMD%"=="clean"      goto :cmd_clean
if /i "%CMD%"=="help"       goto :usage
if /i "%CMD%"=="--help"     goto :usage
if /i "%CMD%"=="-h"         goto :usage
echo [fail] unknown command: %CMD%
call :usage
exit /b 2

:usage
echo Usage: run.bat [verify ^| full [--fast] [--jobs N] ^| app ^| export ^| fetch-data [--yes] ^| clean]
echo        run.bat            (no command) = quick verify, then start the twin
exit /b 0

rem ---------------------------------------------------------------------------
rem  Environment bootstrap: create .venv once, install when requirements change
rem ---------------------------------------------------------------------------
:bootstrap
if exist "%PY%" goto :install
echo [run] creating .venv (Python 3.13 preferred, 3.12 as fallback)
py -3.13 -m venv "%VENV%" 2>nul
if exist "%PY%" goto :install
py -3.12 -m venv "%VENV%" 2>nul
if exist "%PY%" goto :install
echo [warn] Python 3.13 / 3.12 not found through the "py" launcher.
echo [warn] Falling back to the default "python". The study is tested on 3.13 and 3.12 only.
python -m venv "%VENV%"
if exist "%PY%" goto :install
echo [fail] could not create .venv. Install Python 3.13 from https://www.python.org/downloads/
exit /b 1

:install
rem Python itself compares the hash of requirements.txt with the stamp left by the last install
rem (exit code 0 = unchanged, so the install is skipped).
"%PY%" -c "import hashlib,pathlib,sys; h=hashlib.sha256(open('requirements.txt','rb').read()).hexdigest(); p=pathlib.Path(sys.argv[1]); sys.exit(0 if p.exists() and p.read_text().strip()==h else 1)" "%STAMP%"
if not errorlevel 1 exit /b 0
echo [run] installing Python packages from requirements.txt (first run, or the file changed)
"%PY%" -m pip install --disable-pip-version-check -r requirements.txt
if errorlevel 1 (
  echo [fail] package installation failed
  exit /b 1
)
"%PY%" -c "import hashlib,pathlib,sys; pathlib.Path(sys.argv[1]).write_text(hashlib.sha256(open('requirements.txt','rb').read()).hexdigest())" "%STAMP%"
echo [ ok] packages installed
exit /b 0

rem ---------------------------------------------------------------------------
rem  Commands
rem ---------------------------------------------------------------------------
:cmd_default
call :bootstrap
if errorlevel 1 exit /b 1
"%PY%" -m eigenmap.cli verify
if errorlevel 1 (
  echo [fail] verification failed ^(see the table above^) - the twin was not started.
  echo [fail] Start it anyway with:  run.bat app
  exit /b 1
)
echo [ ok] verification passed
goto :run_app

:cmd_verify
call :bootstrap
if errorlevel 1 exit /b 1
"%PY%" -m eigenmap.cli verify %ARGS%
exit /b %errorlevel%

:cmd_full
call :bootstrap
if errorlevel 1 exit /b 1
"%PY%" -m eigenmap.cli full %ARGS%
exit /b %errorlevel%

:cmd_export
call :bootstrap
if errorlevel 1 exit /b 1
"%PY%" -m eigenmap.cli export %ARGS%
exit /b %errorlevel%

:cmd_app
call :bootstrap
if errorlevel 1 exit /b 1
goto :run_app

:run_app
if not exist "twin\app.py" (
  echo [warn] the 3D twin ^(twin\app.py^) is not written yet - nothing to launch.
  exit /b 0
)
echo [run] starting the twin at http://127.0.0.1:%PORT%/   (stop it with Ctrl-C)
rem Open the browser about 4 seconds later, in the background (ping is used as a timer).
if not defined EIGENMAP_NO_BROWSER start "" /b cmd /c "ping -n 5 127.0.0.1 >nul & start http://127.0.0.1:%PORT%/"
"%PY%" -m twin.app --port %PORT%
exit /b %errorlevel%

:cmd_clean
rem Needs no packages, so do not build .venv just to clean.
if exist "%PY%" (
  "%PY%" -m eigenmap.cli clean
) else (
  python -m eigenmap.cli clean
)
exit /b %errorlevel%

:cmd_fetch
rem OPTIONAL. The study does not need this data; every failure here is non-fatal (exit 0).
set "URL=https://phm-datasets.s3.amazonaws.com/NASA/17.+Turbofan+Engine+Degradation+Simulation+Data+Set+2.zip"
set "PAGE=https://www.nasa.gov/intelligent-systems-division/discovery-and-systems-health/pcoe/pcoe-data-set-repository/"
set "TARGET=data\N-CMAPSS_Turbofan_Engine_Degradation_Simulation_Data_Set_2.zip"
echo [run] NASA N-CMAPSS: "Turbofan Engine Degradation Simulation Data Set-2" (NASA Prognostics Data Repository)
echo [warn] the download address below could NOT be verified when this script was written:
echo          %URL%
echo        If it fails, download by hand from the repository page (data set no. 17):
echo          %PAGE%
echo        and put the zip file into the folder: data\
echo [warn] the file is LARGE (several GB). The study and the twin do NOT need it.
if /i "%~2"=="--yes" goto :fetch_go
set "ANSWER="
set /p "ANSWER=Download now? [y/N] "
if /i not "%ANSWER%"=="y" (
  echo [run] download skipped.
  exit /b 0
)
:fetch_go
if not exist data mkdir data
curl.exe -fL -C - --retry 2 -o "%TARGET%.part" "%URL%"
if errorlevel 1 (
  echo [warn] download failed ^(this is not an error for the study^). Manual download page:
  echo          %PAGE%
  exit /b 0
)
move /y "%TARGET%.part" "%TARGET%" >nul
echo [ ok] saved %TARGET%
exit /b 0
