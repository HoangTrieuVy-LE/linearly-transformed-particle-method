!---------------------------------------------------------------------------  
!> @author 
!> HoangTrieuVy-LE.
!
! DESCRIPTION: 
!> Name: data_launch
!> Store data set-up subroutine 
!===========================================================================================


module data_launch
	USE PACKMPI
	USE calculsfor_rec_modf90  !Update also these variables 
	use calculsfor_ini_modf90  !Update also these variables 

	IMPLICIT NONE
	

    INTEGER                                        :: time_scheme
	DOUBLE PRECISION                               :: time_step, T_start, T_end
	
	! Simution region limits
	DOUBLE PRECISION , ALLOCATABLE, DIMENSION(:,:) :: mesh
	
	! Particles Coordinates array
	DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)  :: Xparticle_read 

	! Particles deformation matrix
	DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)  :: D_read
	
	! Particles mass array
	DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)    :: mass_read 
	
	INTEGER, DIMENSION(MPI_STATUS_SIZE)            :: status
	
	! Number of particles simulated
	INTEGER                                        :: number_of_particles


CONTAINS
SUBROUTINE setup_variables
	!> Input: 	'trunk/fortran_srcs/size.data', size.data made from initialization fonction in "code_LTP_2D_diffusion.py"  
	!> Output:  Assert values for: nxg, nyg, hx, time_step, T_start, T_end, time_scheme, degre_spline, radius_spline, phi_radial
	
	OPEN(unit = 1, FILE = 'trunk/fortran_srcs/size.data',status ='old')
	READ(1,*) nxg              	! Nx particles 
	READ(1,*) nyg				! Ny particles
	READ(1,*) hx				! hx remap
	READ(1,*) hy				! hy remap
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
	!> Allocation for:  mass_read, Xparticle_read, D_read, mesh
	number_of_particles = nxg*nyg
	ALLOCATE(mass_read(number_of_particles))
	
	ALLOCATE(Xparticle_read(2,number_of_particles))
	
	ALLOCATE(D_read(4,number_of_particles))
	
	ALLOCATE(mesh(2,2))

END SUBROUTINE setup_particle_information_matrix

SUBROUTINE set_Mesh
	!> Input: 'trunk/fortran_srcs/mesh.bin', made from initialization fonction in "code_LTP_2D_diffusion.py" 
	!> Output: Assert value for mesh array
	OPEN(unit = 5, FILE = 'trunk/fortran_srcs/mesh.bin',FORM="UNFORMATTED",&
	STATUS="UNKNOWN", ACTION="READ", ACCESS='STREAM')
	READ(5) mesh
!	print*,'mesh',mesh(1,1),mesh(1,2),mesh(2,1),mesh(2,2)
	CLOSE(5)
END SUBROUTINE set_Mesh

SUBROUTINE mass_initiation
	!> Input: 'trunk/fortran_srcs/matrixM4fortran.bin', made from initialization fonction in "code_LTP_2D_diffusion.py" 
	!> Output: Assert value for mass matrix
	OPEN(unit = 4, FILE = 'trunk/fortran_srcs/matrixM4fortran.bin',FORM="UNFORMATTED",&
	STATUS="UNKNOWN", ACTION="READ", ACCESS='STREAM')
	READ(4) mass_read
	mpart = mass_read
	CLOSE(4)
	
END SUBROUTINE mass_initiation
	

SUBROUTINE particle_coordinates_initiation
	!> Input: 'trunk/fortran_srcs/coords4fortran.bin', made from initialization fonction in "code_LTP_2D_diffusion.py" 
	!> Output: Assert value for Xparticle_read array
	OPEN(unit = 2, FILE = 'trunk/fortran_srcs/coords4fortran.bin',FORM="UNFORMATTED", &
 STATUS="UNKNOWN",POSITION="REWIND", ACTION="READWRITE", ACCESS='STREAM')
	READ(2) Xparticle_read
END SUBROUTINE particle_coordinates_initiation	


SUBROUTINE deformation_matrix_initiation
	!> Input: 'trunk/fortran_srcs/deformmatrix4fortran.bin', made from initialization fonction in "code_LTP_2D_diffusion.py" 
	!> Output: Assert value for D_read array
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
