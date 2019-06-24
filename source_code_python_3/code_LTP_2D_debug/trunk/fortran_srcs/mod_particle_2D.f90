!---------------------------------------------------------------------------  
!> @author 
!> HoangTrieuVy-LE.
!
! DESCRIPTION: 
!> Particle model
!--------------------------------------------------------------------------- 

module PARTICLE_2D_modf90

!import module, TODO decide which module would be necessary
use calculsfor_rec_modf90
use calculsfor_var_modf90
use calculsfor_ini_modf90


implicit none 
	! DECLARATIONS
	! TODO  



type PARTTYPE

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

end type PARTTYPE

contains
!- Particle structure
!type(PARTTYPE), dimension(:,:), allocatable :: PART



	subroutine MPI_PART_TYPE
	

	end subroutine MPI_PART_TYPE

end module PARTICLE_2D_modf90




