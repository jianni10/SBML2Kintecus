@echo off
rem runkintpack.bat - unpack a .kintpack and run Kintecus on it (the .kintpack is never changed).
rem   runkintpack.bat PACK.kintpack [-k "switches"] [-p Kintecus.exe] [-i] [--set sheet.field=value ...]
rem   runkintpack.bat --help        for all options.   Needs Python 3.8+ (no extra packages).
setlocal
set "PY="
where py >nul 2>nul && set "PY=py -3"
if "%PY%"=="" where python >nul 2>nul && set "PY=python"
if "%PY%"=="" (
    echo Python 3 was not found. Install it from https://www.python.org/downloads/ and run this again.
    exit /b 2
)
%PY% "%~dp0runkintpack.py" %*
exit /b %ERRORLEVEL%
