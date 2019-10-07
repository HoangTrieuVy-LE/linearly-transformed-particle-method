1. Down MS-MPI 
2. Create "lib" folder in C:\
3. Copy "MPI" from Mircosoft SDKs to "lib"
4. Copy msmpi.dll from system32 to lib
5. Copy msmpi.lib from MPI\Lib\x64
5. https://abhilashreddy.com/writing/3/mpi_instructions.html
6. Create a "new" "lib" in MPI, copy libmsmpi.a in



7. compile on C/C++ 
IMPORTANTTTTTTTTTTTTT: include <'C:\lib\MPI\Include\mpi.h'>

		gcc -o out.exe my_mpi_program.c -IC:\lib\mpi\include -LC:\lib\mpi\lib -lmsmpi -fno-range-check
		mpiexec -n 4 out.exe 

8. compile on gfortran
		gfortran -o out.exe my_mpi_program.c -IC:\lib\mpi\include -LC:\lib\mpi\lib -lmsmpi -fno-range-check
		mpiexec -n 4 out.exe	
IMPORTANTTTTTTTTTTTTT: include 'C:\lib\MPI\Include\mpif.h'
			include 'C:\lib\MPI\Include\mpif.h'
			include 'C:\lib\MPI\Include\mpif.h'
			include 'C:\lib\MPI\Include\mpif.h'


9999. http://www.math.ucla.edu/~wotaoyin/windows_coding.html