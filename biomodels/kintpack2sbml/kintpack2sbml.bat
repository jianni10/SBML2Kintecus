@echo off
rem kintpack2sbml.bat - export a Kintecus model (.kintpack or deck folder) to SBML / SED-ML / Kintecus XML.
rem   kintpack2sbml.bat INPUT [options]      kintpack2sbml.bat --help for all options
rem   Needs Python 3.8+ and, for SBML output:  pip install python-libsbml
setlocal
set "PY="
where py >nul 2>nul && set "PY=py -3"
if "%PY%"=="" where python >nul 2>nul && set "PY=python"
if "%PY%"=="" (
    echo Python 3 was not found. Install it from https://www.python.org/downloads/ and run this again.
    exit /b 2
)
%PY% "%~dp0kintpack2sbml.py" %*
exit /b %ERRORLEVEL%
