!---------------------------------------------------------------------------  
!> @author 
!> HoangTrieuVy-LE.
!
! DESCRIPTION: 
!> MPI 2D structures 
!--------------------------------------------------------------------------- 


module MPI_2D_structures_modf90

	! Import modules, TODO decide which module would be necessary
	use calculsfor_rec_modf90
	use calculsfor_var_modf90
	use calculsfor_ini_modf90


implicit none

! DECALARATIONS
! TODO


contains


	subroutine overlap_criterion()
	! Check a particle if it satifies overlap_criterion
	! Retrieve the eigenvalues and eigenvectors to determine the direction movement of particles, 	then 
	! decide it leaves the block or not
	! TODO
	end subroutine overlap_criterion


	subroutine particle_numerotation()
	! Numerize particle inside block
	end subroutine



end module MPI_2D_structures_modf90
