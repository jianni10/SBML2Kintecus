@echo off
rem ======================================================================================
rem  RUN_SLOW_JOBS.bat - run the packs that did not finish within the time limits used for
rem  the BioModels batch (listed in SLOW_JOBS.csv), with NO time limit, each with its own
rem  recommended integrator, and compare each result with its SBML reference.
rem
rem  Success = Kintecus output under 4%% of the SBML (libRoadRunner) result for every
rem  compared species (max |difference| / species peak).
rem
rem  Results:  slow_jobs_report.txt and slow_jobs_report.csv (one row per pack). Every run
rem  folder is kept in kintcheck_work\<ID>\ (CONC.TXT etc.).
rem
rem  Usage:  double-click, or from a command prompt:
rem     RUN_SLOW_JOBS.bat                               all jobs, fastest first
rem     RUN_SLOW_JOBS.bat --only BIOMD0000000017        one job (switches still from SLOW_JOBS.csv)
rem     RUN_SLOW_JOBS.bat --only BIOMD0000000017 --switches "-ig:mass -INT:5 -obeymaxint"
rem     RUN_SLOW_JOBS.bat --cap 86400                   give up on any job after 24 h
rem     RUN_SLOW_JOBS.bat --stall 3600                  give up if no output row for 1 h
rem  Progress (simulated time, % done, minutes to go) is printed every 60 s.
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
%PY% "%~dp0kintcheck.py" --kintecus "%KINTECUS_EXE%" --packs "%~dp0." --refs "%~dp0references" --jobs "%~dp0SLOW_JOBS.csv" --cap 0 --stall 0 --progress 60 --keep --name slow_jobs_report %*
echo.
echo Done. See slow_jobs_report.txt and slow_jobs_report.csv in this folder.
pause
