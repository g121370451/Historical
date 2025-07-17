@echo off
rem === 参数说明 ===
rem %1 = dataset（如 enron-email）
rem %2 = k
rem %3 = t
rem %4 = max
rem %5 = min

setlocal

set BIN_DIR=%~dp0.\bin
set data_DIR=%~dp0.\data
set EXEC=%BIN_DIR%\experiment_program.exe

if "%1"=="" (
    echo [ERROR] Missing parameter: dataset
    goto usage
)
if "%2"=="" (
    echo [ERROR] Missing parameter: k
    goto usage
)
if "%3"=="" (
    echo [ERROR] Missing parameter: t
    goto usage
)
if "%4"=="" (
    echo [ERROR] Missing parameter: max
    goto usage
)
if "%5"=="" (
    echo [ERROR] Missing parameter: min
    goto usage
)

set dataset=%1
set k=%2
set t=%3
set max=%4
set min=%5

if /I "%dataset%"=="enron-email" (
    set file=%data_DIR%/enron-email/raw/Email-Enron.txt
    set out=%data_DIR%/enron-email/processed/
) else (
    echo [ERROR] unknown dataset: %dataset%
    goto usage
)

echo.
echo [INFO] Running label generation...
echo %EXEC% generate-label -t %t% -f %file% -p %out% -k %k% -max %max% -min %min%
"%EXEC%" generate-label -t %t% -f %file% -p %out% -k %k% -max %max% -min %min%

echo.
echo [INFO] Plotting degree distribution with Python...
python "%BIN_DIR%\plot_degree.py" "%out%"

goto :eof

:usage
echo.
echo usage: %~nx0 [dataset] [k-limit] [thread] [max-weight] [min-weight]
echo e.g : %~nx0 enron-email 3 8 100 1
echo.
echo dataset supported:
echo    - enron-email
echo    - ....
