@ECHO OFF
IF "%3"=="" ECHO Please enter the subdirectory to be processed and parameter file (new M, new T) & ECHO Example: .\run.bat test Benders 40 60 & GOTO:EOF

set BINDIR=..\VS2022\x64\Release
set DATADIR=..\..\Data\%1
set PARAMFILE=Parameters\%2.txt
set RESULTDIR=Results\%1_%2_%3_%4
set SUMMARYFILE=market_selection_%1_%2_%3_%4.csv
set EXECUTABLE=MIP.exe

set NEW_M=%3
set NEW_T=%4

IF EXIST %RESULTDIR% rmdir /S /Q %RESULTDIR%
md %RESULTDIR%

ECHO Name,^
Parameters,^
M,^
T,^
Lambda,^
IsPareto,^
Heuristic_CPU,^
Heuristic_LB,^
Heuristic_Revenue,^
Heuristic_Cost,^
Heuristic_nSelected,^
MIP_CPU,^
MIP_LB,^
MIP_UB,^
MIP_Revenue,^
MIP_Cost,^
MIP_nSelected,^
AutoBenders_CPU,^
AutoBenders_LB,^
AutoBenders_UB,^
AutoBenders_Revenue,^
AutoBenders_Cost,^
AutoBenders_nSelected,^
AutoBenders_CutCount,^
Decomposition_CPU,^
Decomposition_LB,^
Decomposition_UB,^
Decomposition_Revenue,^
Decomposition_Cost,^
Decomposition_nFixed,^
Decomposition_nSelected,^
Decomposition_CallbackCPU,^
Decomposition_CallCount,^
Decomposition_CutCount,^
 > %RESULTDIR%\%SUMMARYFILE%

FOR %%f IN (%DATADIR%\*.txt) DO (
	ECHO Started solving %%f with %PARAMFILE%
	REM ECHO %BINDIR%\%EXECUTABLE% %%f %RESULTDIR%\%%~nxf %RESULTDIR%\%SUMMARYFILE% %PARAMFILE% %NEW_M% %NEW_T% 
	%BINDIR%\%EXECUTABLE% %%f %RESULTDIR%\%%~nxf %RESULTDIR%\%SUMMARYFILE% %PARAMFILE% %NEW_M% %NEW_T% > %RESULTDIR%\%%~nxf.log & IF ERRORLEVEL 1 GOTO:EOF
)