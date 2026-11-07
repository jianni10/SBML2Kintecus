@echo off
rem ======================================================================================
rem  RUN_REGRESSION.bat - after changing Kintecus, re-run every pack that PASSED when this
rem  release was built (folder 1_success_under_4pct) and list the ones that no longer pass.
rem
rem  Results:  regression_report_FAILED.csv  only the regressions (one row per pack):
rem                FAILED = no longer under 4% of the SBML (or stalled / too slow / FATAL)
rem                WORSE  = still under 4% but the error grew more than 10x (and above 0.1%)
rem            regression_report.txt         summary + the regression table
rem            regression_report.csv         every pack that was run
rem            kintcheck_work\<ID>\          run folders of the failures (CONC.TXT etc.)
rem
rem  Usage:  double-click, or from a command prompt:
rem     RUN_REGRESSION.bat                            all ~890 passing packs (default limits:
rem                                                   120 s without a new row, 1800 s per pack)
rem     RUN_REGRESSION.bat --only BIOMD0000000535     only the listed model IDs
rem     RUN_REGRESSION.bat --cap 600                  shorter limit per pack
rem  Any options are passed to kintcheck.py (python kintcheck.py --help lists them all).
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
%PY% "%~dp0kintcheck.py" --kintecus "%KINTECUS_EXE%" --packs "%~dp0." --refs "%~dp0references" --folder 1_success_under_4pct --stall 120 --cap 1800 --name regression_report %*
echo.
echo Done. Regressions: regression_report_FAILED.csv   (summary: regression_report.txt)
pause
