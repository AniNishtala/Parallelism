Code must require the use of OpenMP. If you are on windows, you must add  OpenMP to the compiler + linker. 

On linux, use sudo apt update. Then do sudo apt install gcc g++ libomp-dev. 


To compile the cpp file. Run the command, g++ -fopenmp –std:c++11 filename.cpp –o filename_executable

Then execute the executable generated from the previous command with the following. ./filename_executable


Command Breakdown is as follows:

g++ is the compiler, -fopenmp is the flag to be added, -std:c++11 specifies the version of C++ to use, filename.cpp is the cpp file to compile, -o is the parameter for output, filename_executable is the executable file generated from the command. 

