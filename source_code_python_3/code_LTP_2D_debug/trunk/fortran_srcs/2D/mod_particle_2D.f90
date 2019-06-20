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


type PARTYPE
!Stored particle informations
!1. Velocity of k-th particle U[1,k],U[2,k]
!2. Deformation Matrix of k-th particle Dk: TODO find eigenvalues, eigenvector???? 
end type PARTYPE

contains

	subroutine MPI_PART_TYPE
	
	end subroutine MPI_PART_TYPE

end module PARTICLE_2D_modf90
