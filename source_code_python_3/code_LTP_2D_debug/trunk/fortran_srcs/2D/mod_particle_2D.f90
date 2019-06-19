!---------------------------------------------------------------------------  
!> @author 
!> HoangTrieuVy-LE.
!
! DESCRIPTION: 
!> Particle model
!--------------------------------------------------------------------------- 

module PARTICLE_2D

!import module, TODO decide which module is necessary
use calculsfor_rec_modf90
use calculsfor_var_modf90
use calculsfor_ini_modf90

implicit none 
	! DECLARATIONS
	! TODO  
type PARTYPE
!Stored particle informations
!1. Velocity U
!2. Deformation Matrix D: TODO find eigenvalues, eigenvector
 
end type PARTYPE

contains

	subroutine MPI_PART_TYPE
	
	end subroutine MPI_PART_TYPE

end module PARTICLE_2D
