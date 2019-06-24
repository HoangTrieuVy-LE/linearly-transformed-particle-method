PROGRAM launch
	IMPLICIT NONE

	include	 'mpif.h'

	INTEGER, DIMENSION(MPI_STATUS_SIZE)           :: status
	INTEGER :: Nx,Ny
	DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: Xparticle 
	DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: D_old
	DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)   :: M_old


	OPEN(unit = 1, FILE = 'trunk/fortran_srcs/size.data',status ='old')
	READ(1,*) Nx
	READ(1,*) Ny
	CLOSE(1)
	print*, Nx
	print*, Ny
	
	ALLOCATE(Xparticle(2,Nx*Ny))
	
	OPEN(unit = 2, FILE = 'trunk/fortran_srcs/coords4fortran.bin',FORM="UNFORMATTED", STATUS="UNKNOWN", ACTION="READ", ACCESS='STREAM')
	READ(2) Xparticle
	CLOSE(2)
	


	
	ALLOCATE(M_old(Nx*Ny))
	OPEN(unit = 4, FILE ='trunk/fortran_srcs/matrixM4fortran.bin',FORM="UNFORMATTED", STATUS="UNKNOWN", ACTION="READ", ACCESS='STREAM')
	READ(4) M_old
	CLOSE(4)
	
	
	
	ALLOCATE(D_old(4,Nx*Ny))
	OPEN(unit= 3, FILE='trunk/fortran_srcs/deformmatrix4fortran.bin',FORM="UNFORMATTED",STATUS="UNKNOWN",ACTION="READ",ACCESS='STREAM')
	READ(3) D_old
	CLOSE(4)

	

	DEALLOCATE(M_old)
	DEALLOCATE(D_old)
	DEALLOCATE(Xparticle)
END PROGRAM launch
