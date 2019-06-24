!---------------------------------------------------------------------------  
!> @author 
!> HoangTrieuVy-LE.
!
! DESCRIPTION: 
!> Tri_caiser_method ? Renamne TODO
!--------------------------------------------------------------------------- 

module tri_casier_modf90


use MPI_2D_structures_modf90
use PARTICLE_2D_modf90

implicit none

! DECLARATIONS TODO

! setup
! call environnement_initialisation
! call mpi_2d_spatial_decomposition
! call adaptation_block_processor

contains
	
	! For each processor, define "table_particle_inside" and "table_particle_overlap"
	
	! table_particle_inside = TODO
	
	! REALLY NEED TO CLARIFY
	! type(PARTTYPE), dimension (:), allocatable :: table_particle_inside
	! table_particle_overlap(NPmax) 
	! table_particle_overlap = TODO
	! type(PARTTYPE), dimension (:,:), allocatable :: table_particle_overlap
	! table_particle_overlap(NPmax, nb_neighbours)
	! Which is a table would show overlap particlesin each neighbour.
	
	subroutine block_loop
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
	
	end subroutine block_loop

end module tri_casier_modf90
