@echo off
rem ======================================================================================
rem  RUN_ALL_CHECKS.bat - run every .kintpack in this folder through Kintecus.exe and
rem  compare each result with the SBML reference output in the "references" folder.
rem
rem  Success = Kintecus output under 4%% of the SBML (libRoadRunner) result for every
rem  compared species (max |difference| / species peak).
rem
rem  Results:  kintcheck_report.txt (summary) and kintcheck_report.csv (one row per pack)
rem
rem  Usage:  double-click, or from a command prompt:
rem     RUN_ALL_CHECKS.bat                          all packs
rem     RUN_ALL_CHECKS.bat --folder 1_kintecus      only packs in folders matching the text
rem     RUN_ALL_CHECKS.bat --only BIOMD0000000010   only the listed model IDs
rem     RUN_ALL_CHECKS.bat --stall 120 --cap 1800   more patience per pack (seconds)
rem  Any options are passed to kintcheck.py (python kintcheck.py --help lists them all).
rem
rem  Needs Python 3.8 or newer (python.org); no extra packages.
rem ======================================================================================
setlocal
cd /d "%~dp0"

rem ---- where is Kintecus.exe?  Edit this line, or set KINTECUS_EXE before running.
if "%KINTECUS_EXE%"=="" set "KINTECUS_EXE=C:\Kintecus\Kintecus.exe"
if not exist "%KINTECUS_EXE%" if exist "%~dp0Kintecus.exe" set "KINTECUS_EXE=%~dp0Kintecus.exe"
if not exist "%KINTECUS_EXE%" (
    echo Kintecus.exe was not found at "%KINTECUS_EXE%".
    set /p "KINTECUS_EXE=Full path to Kintecus.exe: "
)
if not exist "%KINTECUS_EXE%" (
    echo Still not found - stopping.
    pause
    exit /b 1
)

rem ---- find Python (the py launcher first, then python on PATH)
set "PY="
where py >nul 2>nul && set "PY=py -3"
if "%PY%"=="" where python >nul 2>nul && set "PY=python"
if "%PY%"=="" (
    echo Python 3 was not found. Install it from https://www.python.org/downloads/ and run this again.
    pause
    exit /b 1
)

echo Kintecus: %KINTECUS_EXE%
echo Python:   %PY%
echo.
%PY% "%~dp0kintcheck.py" --kintecus "%KINTECUS_EXE%" --packs "%~dp0." --refs "%~dp0references" %*
echo.
echo Done. See kintcheck_report.txt and kintcheck_report.csv in this folder.
pause
