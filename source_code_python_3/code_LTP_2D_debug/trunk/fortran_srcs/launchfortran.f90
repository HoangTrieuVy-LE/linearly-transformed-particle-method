module data_launch
	USE PACKMPI
	USE calculsfor_rec_modf90
	use calculsfor_ini_modf90

	IMPLICIT NONE
    INTEGER                                        :: time_scheme
	DOUBLE PRECISION                               :: dt, Tini, Tmaximum
	DOUBLE PRECISION , ALLOCATABLE, DIMENSION(:,:) :: mesh


	DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)  :: Xparticle_read 
	DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)  :: D_read
	DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)    :: mass_read 
	DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)  :: velocity_read 
	
	

	INTEGER, DIMENSION(MPI_STATUS_SIZE)            :: status
	INTEGER                                        :: fh1,fh2,fh3&
													, bytes_in_DOUBLE_PRECISION


	INTEGER                                        :: number_of_particles
CONTAINS



SUBROUTINE set_Nx_Ny
	OPEN(unit = 1, FILE = 'trunk/fortran_srcs/size.data',status ='old')
	READ(1,*) nxg
	READ(1,*) nyg
	READ(1,*) hx
	READ(1,*) hy
	READ(1,*) dt
	READ(1,*) Tini
	READ(1,*) Tmaximum
	READ(1,*) time_scheme
	READ(1,*) degre_spline
	READ(1,*) radius_phi
	READ(1,*) phi_radial

	CLOSE(1)
END SUBROUTINE set_Nx_Ny

SUBROUTINE set_length_data_read

	number_of_particles = nxg*nyg
	ALLOCATE(mass_read(number_of_particles))
	
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

SUBROUTINE mass_initiation
	
	OPEN(unit = 4, FILE = 'trunk/fortran_srcs/matrixM4fortran.bin',FORM="UNFORMATTED",&
	STATUS="UNKNOWN", ACTION="READ", ACCESS='STREAM')
	READ(4) mass_read
	mpart = mass_read
	CLOSE(4)
	
END SUBROUTINE mass_initiation

SUBROUTINE velocity_initiation
	OPEN(unit = 7, FILE = 'trunk/fortran_srcs/velocity4fortran.bin',FORM="UNFORMATTED",&
STATUS="UNKNOWN", ACTION="READ", ACCESS='STREAM')
	READ(7) velocity_read
	CLOSE(7)
END SUBROUTINE velocity_initiation
	

SUBROUTINE particle_coordinates_initiation
	OPEN(unit = 2, FILE = 'trunk/fortran_srcs/coords4fortran.bin',FORM="UNFORMATTED", &
 STATUS="UNKNOWN",POSITION="REWIND", ACTION="READWRITE", ACCESS='STREAM')
	READ(2) Xparticle_read



END SUBROUTINE particle_coordinates_initiation	


SUBROUTINE deformation_matrix_initiation	
	OPEN(unit= 3, FILE='trunk/fortran_srcs/deformmatrix4fortran.bin',FORM="UNFORMATTED",STATUS="UNKNOWN",POSITION="REWIND" & 
,ACTION="READWRITE",ACCESS='STREAM')
	READ(3) D_read

END SUBROUTINE deformation_matrix_initiation


SUBROUTINE dealloc_XMD
	DEALLOCATE(mesh)
	DEALLOCATE(mass_read)
	DEALLOCATE(D_read)
	DEALLOCATE(Xparticle_read)
END SUBROUTINE dealloc_XMD

END module data_launch
