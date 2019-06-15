Open Anaconda Prompt
Run code_LTP_2D:
	1. cd Users\Eway\Downloads\source_code_particles\code_LTP_2D_debug\trunk\python_srcs
	2. python launchRuns.py data/inputs/LTP_fastest test

To compile fortran modules :

1) cd ./trunk/python_srcs/

2) python setup.py build_ext --inplace


	
Run with MPI: 
c	ompile on gfortran
		gfortran -o out.exe my_mpi_program.c -IC:\lib\mpi\include -LC:\lib\mpi\lib -lmsmpi -fno-range-check
		mpiexec -n 4 out.exe	
IMPORTANTTTTTTTTTTTTT: include 'C:\lib\MPI\Include\mpif.h'
			include 'C:\lib\MPI\Include\mpif.h'
			include 'C:\lib\MPI\Include\mpif.h'
			include 'C:\lib\MPI\Include\mpif.h'