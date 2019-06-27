!---------------------------------------------------------------------------  
!> @author 
!> HoangTrieuVy-LE.
!
! DESCRIPTION: 
!> Tri_caiser_method
!--------------------------------------------------------------------------- 

! TODO setup
! call environnement_initialisation
! call mpi_2d_spatial_decomposition
! call adaptation_block_processor

PROGRAM tri_casier_modf90


USE data_launch
USE mpi_2D_structures_modf90
USE mod_particle_2D_modf90
USE mpi_2D_spatialisation_modf90


IMPLICIT NONE
	logical :: overlapcheck ! case test
! DECLARATIONS TODO
	
	
	!-------------------------------------------------------------------!
	!!!                    PARTICLES INITIALISATION                   !!!
	!-------------------------------------------------------------------!
	CALL environnement_initialisation

	CALL set_Nx_Ny
	CALL set_length_data_read
	CALL set_Mesh
	

	!WRITE(*,*) 'mesh:',mesh(:,1)	
	!WRITE(*,*) 'Nx = ', Nx
	!WRITE(*,*) 'Ny = ', Ny
	!WRITE(*,*) 'length_read', length_read
	
	CALL particle_coordinates_initiation
	!WRITE(*,*) Xparticle_read(:,2)
	CALL density_initiation
	!write(*,*) velocity_read(:,1:10)tion
	!WRITE(*,*) density_read(1:10)
	CALL deformation_matrix_initiation
	CALL velocity_initiation

	 
	!-------------------------------------------------------------------!
	!!!                         DOMAIN CREATION                       !!!
	!-------------------------------------------------------------------!
	CALL topology_initialisation
	!write(*,*) 'Lx1: ', mesh(1,1)
	!write(*,*) 'Lx2: ', mesh(1,2)
	!write(*,*) 'Ly1: ', mesh(2,1)
	!write(*,*) 'Ly2: ', mesh(2,2)
	call set_block_grid(mesh(1,1),mesh(1,2),mesh(2,1),mesh(2,2))

	call neighbour_blocks
	call neighbour_counter
	!print*,'rank:',rank,'nb_neighbours',nb_neighbours	
	!if (rank==3) then 	
!		call neighbour_blocks
!		write(*,*)'rank_left',rank_left
!		write(*,*)'rank_right',rank_right 
!		write(*,*)'rank_up',rank_up
!		write(*,*)'rank_down',rank_down
!		write(*,*)'rank_up_left',rank_up_left 
!		write(*,*)'rank_up_right',rank_up_right
!		write(*,*)'rank_down_left',rank_down_left 
!		write(*,*)'rank_down_right',rank_down_right
	!end if

	!call neighbour_blocks


	!-------------------------------------------------------------------!
	!!!  PARTICLES ATTRIBUTION  &  OVERLAP PARTICLES LISTS CREATION   !!!
	!-------------------------------------------------------------------!
	call initiation_table
	call parttype_convert
	! write(*,*) , 'density : ', ALL_PARTICLES(1)%Rhop
	!print*, 'rank:',rank,'star_x', start_x, 'end_x', end_x, 'start_y', start_y,'end_y',end_y



! For each processor, define "table_particle_inside" and "table_particle_overlap"
	call particle_distribution
	write(*,*)'rank:',rank, 'COUNTER_inside', COUNTER_inside
	write(*,*)'rank:',rank, 'COUNTER_overlap', COUNTER_overlap

!	call overlap_criterion(ALL_PARTICLES(1),overlapcheck)
!	print*,overlapcheck
	
	!-------------------------------------------------------------------!
	!!!          FOR EVERY SUB-DOMAIN - INSIDE A PROCESS              !!!
	!-------------------------------------------------------------------!
	!-------------------------------------------------------------------!
	!!!LOOP ON PARTICLES INSIDE BLOCK AND THE OTHERS IN OVERLAP TABLES!!!
	!-------------------------------------------------------------------! 
	
	! TODO BLOCK LOOP
	call block_loop
	
	!-------------------------------------------------------------------!
	!!!                            UPDATE                             !!!
	!-------------------------------------------------------------------!
	
	! Displacement
	call update_displacement
	! Update all new position
	
	
	! Update table
	! TODO
	
	CALL dealloc_XMD
	CALL environnement_finalization
END PROGRAM tri_casier_modf90
