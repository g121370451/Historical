@echo off
rem === 参数说明 ===
rem %1 = dataset（如 enron-email）
rem %2 = k
rem %3 = t
rem %4 = max
rem %5 = min
rem %6 = m 迭代次数
rem %7 = c 变化次数
rem %8 = strategy 变化策略
rem %9 = (可选) --check

setlocal

set BIN_DIR=%~dp0.\bin
set DATA_DIR=%~dp0.\data
set EXEC=%BIN_DIR%\experiment_program.exe

rem ==== 参数校验 ====
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
if "%6"=="" (
    echo [ERROR] Missing parameter: m
    goto usage
)
if "%7"=="" (
    echo [ERROR] Missing parameter: c
    goto usage
)
if "%8"=="" (
    echo [ERROR] Missing parameter: strategy
    goto usage
)

set dataset=%1
set k=%2
set t=%3
set max=%4
set min=%5
set m=%6
set c=%7
set strategy=%8
set check=%9

rem ==== 支持的策略列表 ====
set STRATEGIES=high_high_increase high_high_decrease high_high_mixed low_low_increase low_low_decrease low_low_mixed high_low_increase high_low_decrease high_low_mixed

rem ==== 检查 strategy 合法性 ====
set found=0
for %%s in (%STRATEGIES%) do (
    if /I "%strategy%"=="%%s" (
        set found=1
    )
)
if "%found%"=="0" (
    echo [ERROR] Invalid strategy: %strategy%
    echo Valid options: %STRATEGIES%
    goto usage
)

if /I "%check%"=="check" (
    set check=--%check%
)

rem ==== 根据数据集设置路径 ====
if /I "%dataset%"=="enron-email" (
    set file=%DATA_DIR%/enron-email/processed/
    set out=%DATA_DIR%/enron-email/processed/
) else (
    echo [ERROR] Unknown dataset: %dataset%
    goto usage
)

rem ==== 执行程序 ====
echo.
echo [INFO] Running label maintain...
echo %EXEC% maintain-label -t %t% -f %file% -p %out% -k %k% -c "%c%" -m "%m%" -max %max% -min %min% --strategy %strategy% %check% -n "%dataset%"
%EXEC% maintain-label -t %t% -f %file% -p %out% -k %k% -c %c% -m %m% -max %max% -min %min% --strategy %strategy% %check% -n "%dataset%"
goto :eof

:usage
echo.
echo usage: %~nx0 [dataset] [k-limit] [thread] [max-weight] [min-weight] [iteration_count] [edge_change_count] [strategy] [check]
echo e.g : %~nx0 enron-email 3 8 100 1 10 500 high_high_increase check
echo.
echo dataset supported:
echo    - enron-email
echo strategy supported:
echo    - high_high_increase
echo    - high_high_decrease
echo    - high_high_mixed
echo    - mid_mid_increase
echo    - mid_mid_decrease
echo    - mid_mid_mixed
echo    - low_low_increase
echo    - low_low_decrease
echo    - low_low_mixed
