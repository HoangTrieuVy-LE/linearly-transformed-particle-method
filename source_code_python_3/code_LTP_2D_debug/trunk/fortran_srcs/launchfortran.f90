module data_launch
	USE PACKMPI
	USE calculsfor_rec_modf90  !Update also these variables 
	use calculsfor_ini_modf90  !Update also these variables 

	IMPLICIT NONE
    INTEGER                                        :: time_scheme
	DOUBLE PRECISION                               :: time_step, T_start, T_end
	DOUBLE PRECISION , ALLOCATABLE, DIMENSION(:,:) :: mesh
	DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)  :: Xparticle_read 
	DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)  :: D_read
	DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)    :: mass_read 
	
	INTEGER, DIMENSION(MPI_STATUS_SIZE)            :: status
	INTEGER                                        :: number_of_particles


CONTAINS
SUBROUTINE setup_variables
	OPEN(unit = 1, FILE = 'trunk/fortran_srcs/size.data',status ='old')
	READ(1,*) nxg
	READ(1,*) nyg
	READ(1,*) hx
	READ(1,*) hy
	READ(1,*) time_step
	READ(1,*) T_start
	READ(1,*) T_end
	READ(1,*) time_scheme
	READ(1,*) degre_spline
	READ(1,*) radius_phi
	READ(1,*) phi_radial

	CLOSE(1)
END SUBROUTINE setup_variables

SUBROUTINE setup_particle_information_matrix

	number_of_particles = nxg*nyg
	ALLOCATE(mass_read(number_of_particles))
	
	ALLOCATE(Xparticle_read(2,number_of_particles))
	
	ALLOCATE(D_read(4,number_of_particles))
	
	ALLOCATE(mesh(2,2))

END SUBROUTINE setup_particle_information_matrix

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
