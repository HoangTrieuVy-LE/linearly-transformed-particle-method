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
USE mpi_2D_tools_modf90
USE mod_particle_2D_modf90
USE mpi_2D_spatialisation_modf90


IMPLICIT NONE

! DECLARATIONS TODO
	
	
	!-------------------------------------------------------------------!
	!!!                    PARTICLES INITIALISATION                   !!!
	!-------------------------------------------------------------------!
	CALL set_Nx_Ny
	
	WRITE(*,*) 'Nx = ', Nx
	WRITE(*,*) 'Ny = ', Ny
	
	CALL particle_coordinates_initiation
	CALL density_initiation
	CALL deformation_matrix_initiation
	
	 
	!-------------------------------------------------------------------!
	!!!                         DOMAIN CREATION                       !!!
	!-------------------------------------------------------------------!
	
	
	 
	!-------------------------------------------------------------------!
	!!!                      PARTICLES ATTRIBUTION                    !!!
	!-------------------------------------------------------------------!
	
	
	
	!-------------------------------------------------------------------!
	!!!                 OVERLAP PARTICLES LISTS CREATION              !!!
	!-------------------------------------------------------------------!
	
	
	! For each processor, define "table_particle_inside" and "table_particle_overlap"
	
	
	! table_particle_inside = TODO
	! REALLY NEED TO CLARIFY
	! type(PARTTYPE), dimension (:), allocatable :: table_particle_inside
	! table_particle_overlap(NPmax) 
	
	
	! table_particle_overlap = TODO
	! type(PARTTYPE), dimension (:,:), allocatable :: table_particle_overlap
	! table_particle_overlap(NPmax, nb_neighbours)
	! Which is a table would show overlap particlesin each neighbour.
	
	
	!-------------------------------------------------------------------!
	!!!          FOR EVERY SUB-DOMAIN - INSIDE A PROCESS              !!!
	!-------------------------------------------------------------------!
	
	

	!-------------------------------------------------------------------!
	!!!LOOP ON PARTICLES INSIDE BLOCK AND THE OTHERS IN OVERLAP TABLES!!!
	!-------------------------------------------------------------------! 
	
	! TODO BLOCK LOOP
	
	!-------------------------------------------------------------------!
	!!!                            UPDATE                             !!!
	!-------------------------------------------------------------------!
	
	! Displacement
	! TODO
	! Update all new position
	
	
	! Update table
	! TODO
	

END PROGRAM tri_casier_modf90
