!---------------------------------------------------------------------------  
!> @author 
!> HoangTrieuVy-LE.
!
! DESCRIPTION: 
!> Tri_caiser_method
!--------------------------------------------------------------------------- 

PROGRAM tri_casier_modf90

USE mpi_2D_structures_modf90
USE mpi_2D_spatialisation_modf90

IMPLICIT NONE

	!-------------------------------------------------------------------!
	!!!                         INITIALISATION                        !!!
	!-------------------------------------------------------------------!
CALL environnement_initialisation

CALL setup_variables
CALL setup_particle_information_matrix
CALL set_Mesh
CALL particle_coordinates_initiation
CALL mass_initiation
CALL deformation_matrix_initiation


	!-------------------------------------------------------------------!
	!!!                    SHAPE FUNCTIONS INITIATION                 !!!
	!-------------------------------------------------------------------!
	!  ALREADY DONE IN PYTHON WITH DATA_CODE_LTP_2D.py line 194

	!-------------------------------------------------------------------!
	!!!                         DOMAIN CREATION                       !!!
	!-------------------------------------------------------------------!
CALL topology_initialisation

CALL set_block_grid(mesh(1,1),mesh(1,2),mesh(2,1),mesh(2,2))

CALL neighbour_blocks
!	if (rank==1) then
!		call neighbour_counter
!		print*,'rank:',rank,'number_of_neighbours',nb_neighbours_actually	 	
!		write(*,*)'rank_left',rank_left
!		write(*,*)'rank_right',rank_right 
!		write(*,*)'rank_up',rank_up
!		write(*,*)'rank_down',rank_down
!		write(*,*)'rank_up_left',rank_up_left 
!		write(*,*)'rank_up_right',rank_up_right
!		write(*,*)'rank_down_left',rank_down_left 
!		write(*,*)'rank_down_right',rank_down_right
!	end if

	!-------------------------------------------------------------------!
	!!!  PARTICLES ATTRIBUTION  &  OVERLAP PARTICLES LISTS CREATION   !!!
	!-------------------------------------------------------------------!
CALL initiation_table
CALL parttype_convert
CALL neighbouring

CALL particle_distribution

	!-------------------------------------------------------------------!
	!!!          FOR EVERY SUB-DOMAIN - INSIDE A PROCESS              !!!
	!-------------------------------------------------------------------!
	!-------------------------------------------------------------------!
	!!!LOOP ON PARTICLES INSIDE BLOCK AND THE OTHERS IN OVERLAP TABLES!!!
	!-------------------------------------------------------------------!
	


DO WHILE(T_start<=T_end)
	
	if (rank==0) then
		print*,'T_start',T_start
		print*,'T_end',T_end
		print*,'Time_step',time_step
		DO i=0,nb_proc-1
			write(*,*)'rank:',i, 'COUNTER_inside', COUNTER_inside
		END DO
	write(*,*)'              				   FJ	','	BJ	','   FK	','	BK	','   FJBK      ','BJFK	  ','BJBK		','FJFK'
		DO i = 0,nb_proc-1
			write(*,*)'rank:',i, 'COUNTER_overlap', COUNTER_overlap(:,i)
		END DO
		DO i = 0,nb_proc-1
			write(*,*)'rank:',i, 'COUNTER_danger', COUNTER_danger(:,i)
		END DO
	end if 
	if (rank==0) then
	DO i = 0,nb_proc-1
			write(*,*)'rank:',i, 'COUNTER_leave', COUNTER_leave(:,i)
		END DO
	end if


	CALL send_overlap_and_danger_particle
	CALL block_loop_on_block_global_table


	!-------------------------------------------------------------------!
	!!!                            UPDATE                             !!!
	!-------------------------------------------------------------------!
	
	! Displacement and Deformation matrix
	CALL update_ALL_particles
	CALL send_all_block_information_to_rank_0

!	
	T_start = T_start + time_step
	
END DO	
	
	!-------------------------------------------------------------------!
	!!!                        UPDATE FILE IN                         !!!
	!-------------------------------------------------------------------!
	! TODO, update mass_initiation, particle_coordinates_initiation, 
	! deformation_matrix_initiation,
	
!	
call update_all_particle_information
		
!	CALL dealloc_XMD
!	CALL dealloc_all_particle_table
CALL environnement_finalization
END PROGRAM tri_casier_modf90
