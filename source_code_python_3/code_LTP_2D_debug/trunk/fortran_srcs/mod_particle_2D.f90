!---------------------------------------------------------------------------  
!> @author 
!> HoangTrieuVy-LE.
!
! DESCRIPTION: 
!> Particle model
!--------------------------------------------------------------------------- 

MODULE mod_particle_2D_modf90

!import module, TODO decide which module would be necessary

IMPLICIT NONE 
	! DECLARATIONS
	! TODO  

INCLUDE 'mpif.h'

INTEGER :: code


TYPE PARTTYPE

!Stored particle informations

!0. Coordinates
real(kind=8) :: Xp
real(kind=8) :: Yp

!1. Velocity of k-th particle U[1,k],U[2,k]
real(kind=8) :: Upx
real(kind=8) :: Upy
 

!2. Deformation Matrix of k-th particle Dk: TODO find eigenvalues, eigenvector????
real(kind=8) :: Dp1 
real(kind=8) :: Dp2
real(kind=8) :: Dp3
real(kind=8) :: Dp4

!3. Poids
real(kind=8) :: Mp

!4. Rho
real(kind=8) :: Rhop

!5. Particle dentity
integer :: ID
integer :: PROC_ID

END TYPE PARTTYPE

!integer :: NINTEGER = 2
integer :: NINTEGER = 0

! NDOUBLE 
integer :: NDOUBLE = 28 + 2500*11

!! MPI TYPE FOR PARTICLE EXCHANGE
integer :: MPI_PARTICLETYPE
integer :: TAG



!- Particle structure
TYPE(PARTTYPE), dimension(:,:), allocatable :: PART


CONTAINS

	SUBROUTINE MPI_PART_TYPE

    INTEGER, DIMENSION(0:1) ::  oldtypes, blockcounts, offsets(0:1)
!    integer(kind=MPI_OFFSET_KIND) :: extent
    INTEGER :: extent

!   Setup description of the NDOUBLE MPI_DOUBLE_PRECISION fields in PARTTYPE  
    offsets(0) = 0 
    oldtypes(0) = MPI_DOUBLE_PRECISION
    blockcounts(0) = NDOUBLE

!  Setup description of the NINTEGER MPI_INTEGER fields in PARTYPE
!  Need to first figure offset by getting size of MPI_DOUBLE_PRECISION 
    CALL MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION, extent, code) 
    offsets(1) = NDOUBLE * extent 
    oldtypes(1) = MPI_INTEGER 
    blockcounts(1) = NINTEGER 

!  Now define structured type and commit it  
    CALL MPI_TYPE_STRUCT(2, blockcounts, offsets, oldtypes, MPI_PARTICLETYPE, code) 
    CALL MPI_TYPE_COMMIT(MPI_PARTICLETYPE, code) 

	END SUBROUTINE MPI_PART_TYPE

END MODULE mod_particle_2D_modf90




