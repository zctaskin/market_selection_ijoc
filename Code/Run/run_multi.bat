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
nPareto,^
nDiscovered,^
nImproved,^
nEliminated,^
FrontierCalculator_CPU,^
FrontierCalculator_Initial,^
FrontierCalculator_Pareto,^
MIP_CPU,^
MIP_Iterations,^
MIP_Improved,^
AutoBenders_CPU,^
AutoBenders_Iterations,^
AutoBenders_Improved,^
AutoBenders_CutCount,^
Decomposition_CPU,^
Decomposition_Iterations,^
Decomposition_Improved,^
Decomposition_CallbackCPU,^
Decomposition_CallCount,^
Decomposition_CutCount,^
DecompositionReuse_CPU,^
DecompositionReuse_Iterations,^
DecompositionReuse_Improved,^
DecompositionReuse_CallbackCPU,^
DecompositionReuse_CallCount,^
DecompositionReuse_CutCount,^
 > %RESULTDIR%\%SUMMARYFILE%

FOR %%f IN (%DATADIR%\*.txt) DO (
	ECHO Started solving %%f with %PARAMFILE%
	REM ECHO %BINDIR%\%EXECUTABLE% %%f %RESULTDIR%\%%~nxf %RESULTDIR%\%SUMMARYFILE% %PARAMFILE% %NEW_M% %NEW_T% 
	%BINDIR%\%EXECUTABLE% %%f %RESULTDIR%\%%~nxf %RESULTDIR%\%SUMMARYFILE% %PARAMFILE% %NEW_M% %NEW_T% > %RESULTDIR%\%%~nxf.log & IF ERRORLEVEL 1 GOTO:EOF
)