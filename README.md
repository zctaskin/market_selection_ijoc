# Overview
This repository contains the source code, data files and parameters for the paper ``A Decomposition Algorithm for Single and Multi-Objective Integrated Market Selection and Production Planning" by Wilco van den Heuvel, Semra Ağralı and Z. Caner Taşkın.

# Dependencies
* Microsoft Visual Studio 2022 with C++ language support
* IBM ILOG CPLEX 20.1 Windows x64 version
# File structure
* Data/: Data files. Each problem instance has a corresponding .txt file. File format is given as comments starting with # in each file. 
* Code: Source code 
* Code\VS2022: Visual Studio project and solution files
* Code\Run: Batch files to run experiments for Tables 2 - 6 
* Code\Parameters: Parameter files for computational experiments

# Build and run instructions
* Clone repository to a local folder. Ensure that dependencies are installed in their default locations (if not, adjust Visual Studio project settings to match the installed location and CPLEX version).
* Open Code\VS2022\MIP.sln file in Visual Studio.
* Choose ``Release" build configuration and build the solution.
* Navigate to Code\Run on command line
* Run table2.bat, table3.bat, table4.bat, table56.bat to run computational experiments reported in the paper. 
* Outputs will be generated in Code\Results folder. 
