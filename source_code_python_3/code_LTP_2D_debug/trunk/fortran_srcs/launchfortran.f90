module data_launch
	USE PACKMPI
	IMPLICIT NONE
	
	DOUBLE PRECISION , ALLOCATABLE, DIMENSION(:,:) :: mesh
	INTEGER                                        :: Nx,Ny
	DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)  :: Xparticle_read 
	DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)  :: D_read
	DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)    :: density_read 
	DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)  :: velocity_read 
	
	

	INTEGER, DIMENSION(MPI_STATUS_SIZE)            :: status
	INTEGER                                        :: fh1,fh2,fh3&
													, bytes_in_DOUBLE_PRECISION


	INTEGER                                        :: number_of_particles
CONTAINS



SUBROUTINE set_Nx_Ny
	OPEN(unit = 1, FILE = 'trunk/fortran_srcs/size.data',status ='old')
	READ(1,*) Nx
	READ(1,*) Ny
	CLOSE(1)
END SUBROUTINE set_Nx_Ny


SUBROUTINE set_length_data_read

	number_of_particles = Nx*Ny
	ALLOCATE(density_read(number_of_particles))
	
	ALLOCATE(Xparticle_read(2,number_of_particles))
	
	ALLOCATE(D_read(4,number_of_particles))
	
	ALLOCATE(mesh(2,2))

	ALLOCATE(velocity_read(2,number_of_particles))

END SUBROUTINE set_length_data_read


SUBROUTINE set_Mesh
	OPEN(unit = 5, FILE = 'trunk/fortran_srcs/mesh.bin',FORM="UNFORMATTED",&
	STATUS="UNKNOWN", ACTION="READ", ACCESS='STREAM')
	READ(5) mesh
	CLOSE(5)
END SUBROUTINE set_Mesh



SUBROUTINE density_initiation
	!call MPI_FILE_OPEN(MPI_COMM_WORLD,'trunk/fortran_srcs/matrixM4fortran.bin',&
	! MPI_MODE_RDWR, MPI_INFO_NULL,fh1,code)
	!call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, bytes_in_DOUBLE_PRECISION,code)
	!write(*,*) 'rank = ', rank
	!offset = rank*number_of_particles*bytes_in_DOUBLE_PRECISION

	!call MPI_FILE_READ(fh1,offset,density_read,number_of_particles,bytes_in_DOUBLE_PRECISION,status,code)
	!print*,'SHAPE: ',shape(density_read)
	!call MPI_FILE_CLOSE(fh1,code)
	OPEN(unit = 4, FILE = 'trunk/fortran_srcs/matrixM4fortran.bin',FORM="UNFORMATTED",&
	STATUS="UNKNOWN", ACTION="READ", ACCESS='STREAM')
	READ(4) density_read
	CLOSE(4)
	
END SUBROUTINE density_initiation

SUBROUTINE velocity_initiation
	OPEN(unit = 7, FILE = 'trunk/fortran_srcs/velocity4fortran.bin',FORM="UNFORMATTED",&
STATUS="UNKNOWN", ACTION="READ", ACCESS='STREAM')
	READ(7) velocity_read
	CLOSE(7)
END SUBROUTINE velocity_initiation
	

SUBROUTINE particle_coordinates_initiation
	OPEN(unit = 2, FILE = 'trunk/fortran_srcs/coords4fortran.bin',FORM="UNFORMATTED", STATUS="UNKNOWN", ACTION="READ", ACCESS='STREAM')
	READ(2) Xparticle_read
	CLOSE(2)
	!call MPI_FILE_OPEN(MPI_COMM_WORLD,'trunk/fortran_srcs/coords4fortran.bin',&
	!MPI_MODE_RDWR,MPI_INFO_NULL,fh2,code)
	!call MPI_TYPE_CONTIGUOUS(2,MPI_DOUBLE_PRECISION,type_column,code)
	!call MPI_TYPE_COMMIT(type_column,code)
	!CALL MPI_TYPE_SIZE(type_column, bytes_in_DOUBLE_PRECISION,code)
	!offset = rank*number_of_particles*bytes_in_DOUBLE_PRECISION
	
	!call MPI_FILE_READ_AT(fh2,offset,Xparticle_read,number_of_particles,MPI_DOUBLE_PRECISION,status,code)

	!call MPI_FILE_CLOSE(fh2,code)

END SUBROUTINE particle_coordinates_initiation	


SUBROUTINE deformation_matrix_initiation	
	OPEN(unit= 3, FILE='trunk/fortran_srcs/deformmatrix4fortran.bin',FORM="UNFORMATTED",STATUS="UNKNOWN",ACTION="READ",ACCESS='STREAM')
	READ(3) D_read
	CLOSE(3)
END SUBROUTINE deformation_matrix_initiation


SUBROUTINE dealloc_XMD
	DEALLOCATE(mesh)
	DEALLOCATE(density_read)
	DEALLOCATE(D_read)
	DEALLOCATE(Xparticle_read)
END SUBROUTINE dealloc_XMD

END module data_launch
