The code is to be complied using general fortran code.
Open Cygwin and go to the folder in which you have stored said files.

Serial Code:

gfortran -o run.exe Serial.f90

Serial OpenMP code:

gfortran - fopenmp OpenMP.f90 -o run.exe
export OMP_NUM_THREADS = 2
./ run.exe

Hybrid Code:

mpif90 - fopenmp -o run.exe OpenMP.f90 
export OMP_NUM_THREADS = 2
mpirun -np 4 ./ run.exe

All these codes will generate an output (Output file if run on Artemis) and a Tec plot which can be visualised using Tec Plot.