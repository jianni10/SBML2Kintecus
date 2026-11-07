@echo off
rem ======================================================================================
rem  MAKE_KINTPACK.bat - convert SBML file(s) to Kintecus .kintpack + DAT files.
rem  Drag one or more SBML files onto this .bat, or run:  MAKE_KINTPACK.bat model.xml [options]
rem  Options are passed to sbml2kintpack.py (python sbml2kintpack.py --help), e.g.
rem      MAKE_KINTPACK.bat model.xml --t-end 500 --points 1000
rem  Needs Python 3.8+ and, once:  pip install python-libsbml sympy libroadrunner numpy
rem  If Kintecus.exe is found (path below), long USER expressions are checked in Kintecus.
rem ======================================================================================
setlocal
if "%KINTECUS_EXE%"=="" set "KINTECUS_EXE=C:\Kintecus\Kintecus.exe"
set "PY="
where py >nul 2>nul && set "PY=py -3"
if "%PY%"=="" where python >nul 2>nul && set "PY=python"
if "%PY%"=="" (
    echo Python 3 was not found. Install it from https://www.python.org/downloads/ and run this again.
    pause
    exit /b 1
)
set "KOPT="
if exist "%KINTECUS_EXE%" set KOPT=--kintecus "%KINTECUS_EXE%"
if "%~1"=="" (
    echo Drag SBML files onto this .bat, or run: MAKE_KINTPACK.bat model.xml [options]
    pause
    exit /b 1
)
set "SECOND=%~2"
if "%SECOND:~0,2%"=="--" goto :one
set "EXT=%~x1"
if /i "%EXT%"==".xml" goto :many
if /i "%EXT%"==".sbml" goto :many
goto :one
:many
rem several dropped files: convert each one
for %%F in (%*) do (
    if /i "%%~xF"==".xml" %PY% "%~dp0sbml2kintpack.py" "%%~fF" %KOPT%
    if /i "%%~xF"==".sbml" %PY% "%~dp0sbml2kintpack.py" "%%~fF" %KOPT%
)
goto :end
:one
%PY% "%~dp0sbml2kintpack.py" %* %KOPT%
:end
echo.
pause
