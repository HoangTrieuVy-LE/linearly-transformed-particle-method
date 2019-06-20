!---------------------------------------------------------------------------  
!> @author 
!> HoangTrieuVy-LE.
!
! DESCRIPTION: 
!> Tri_caiser_method ? Renamne TODO
!--------------------------------------------------------------------------- 

module tri_casier_modf90

!import module, TODO decide which module would be necessary

use MPI_2D_spatialisation_modf90
use MPI_2D_structures_modf90
use PARTICLE_2D_modf90

implicit none
s
! DECLARATIONS TODO

contains
	
	! For each processor, define "table_particle_inside" and "table_particle_overlap"
	
	! table_particle_inside = TODO
	! table_particle_overlap = TODO

	subroutine table_particle_inside
	end subroutine table_particle_inside

	subroutine table_particle_overlap
	end subroutine table_particle_overlap


	subroutine main
	!-------------------------------------------------------------------!
	!!!LOOP ON PARTICLES INSIDE BLOCK AND THE OTHERS IN OVERLAP TABLES!!!
	!-------------------------------------------------------------------! 
	! TODO

	! Input: X_casier_k (k is rank of processor),M_casier_k,D_old_casier_k,hx,hy, N_part_casier_k	
	! Output: U_casier_k, Df2py_casier_k
	
	
	! Displacement
	! TODO
	! Update all new position
	
	
	! Update local table
	! TODO
	! Remake step define the two "table_particle_inside" and "table_particle_overlap"	
	
	end subroutine main


end module tri_casier_modf90
