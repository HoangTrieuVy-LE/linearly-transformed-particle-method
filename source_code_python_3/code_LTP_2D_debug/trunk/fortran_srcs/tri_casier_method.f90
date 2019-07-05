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
!	logical :: overlapcheck ! case test
!	double precision, dimension(1,nb_proc) :: matrix
! DECLARATIONS TODO
	
	
	!-------------------------------------------------------------------!
	!!!                    PARTICLES INITIALISATION                   !!!
	!-------------------------------------------------------------------!
	CALL environnement_initialisation

	CALL set_Nx_Ny
	CALL set_length_data_read
	CALL set_Mesh
	
	CALL particle_coordinates_initiation

	CALL mass_initiation

	CALL deformation_matrix_initiation
	CALL velocity_initiation

	!-------------------------------------------------------------------!
	!!!                    SHAPE FUNCTIONS INITIATION                 !!!
	!-------------------------------------------------------------------!
	!  ALREADY DONE IN PYTHON WITH DATA_CODE_LTP_2D.py line 194

	!-------------------------------------------------------------------!
	!!!                         DOMAIN CREATION                       !!!
	!-------------------------------------------------------------------!
	CALL topology_initialisation

	call set_block_grid(mesh(1,1),mesh(1,2),mesh(2,1),mesh(2,2))

	call neighbour_blocks
!	call neighbour_counter
!	print*,'rank:',rank,'nb_neighbours',nb_neighbours	
!	if (rank==1) then 	
!		call neighbour_blocks
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
	call initiation_table
	call parttype_convert
	
	call particle_distribution
	
	if (rank==0) then
	DO i=0,nb_proc-1
		write(*,*)'rank:',i, 'COUNTER_inside', COUNTER_inside
	END DO
	write(*,*)'              				   FJ	','	BJ	','   FK	','	BK	','   FJBK	      ','BJFK	  ','BJBK		','FJFK		'
	DO i = 0,nb_proc-1	
		write(*,*)'rank:',i, 'COUNTER_overlap', COUNTER_overlap(:,i)
	END DO
	DO i = 0,nb_proc-1	
		write(*,*)'rank:',i, 'COUNTER_danger', COUNTER_danger(:,i)
	END DO
	DO i = 0,nb_proc-1	
		write(*,*)'rank:',i, 'COUNTER_leave', COUNTER_leave(:,i)
	END DO

	end if 
	
	call neighbouring
	call send_overlap_and_danger_particle

	
	!-------------------------------------------------------------------!
	!!!          FOR EVERY SUB-DOMAIN - INSIDE A PROCESS              !!!
	!-------------------------------------------------------------------!
	!-------------------------------------------------------------------!
	!!!LOOP ON PARTICLES INSIDE BLOCK AND THE OTHERS IN OVERLAP TABLES!!!
	!-------------------------------------------------------------------!
!	if (rank==0) then
!		print*,'BEFORE', ALL_PARTICLES(1)%Xp
!	end if 

	call block_loop_on_block_global_table
	call update_displacement_on_just_IND_inside_and_IND_overlap

	if (rank==0) then
		print*,'AFTER', ALL_PARTICLES(1)%Xp
	end if 
	
	 

!DO WHILE(Tini<=Tmaximum)
!! TODO
!	!update_overlap_table
!	if(rank==0) then
!	print*, "t= ", Tini
!	end if
!!	if(rank==0) then
!!		write(*,*) 'before', Xparticle_read(1,1:5)
!!	end if
!	call block_loop
!	
!	!-------------------------------------------------------------------!
!	!!!                            UPDATE                             !!!
!	!-------------------------------------------------------------------!
!	
!	! Displacement
!	call update_displacement
!!	call update_deformation_matrix
!	call update_table_block

!!	if(rank==0) then
!!	write(*,*) 'after', Xparticle_read(1,1:5)
!!	end if
!	if (rank==0) then
!		write(*,*)'rank:',rank, 'COUNTER_inside', COUNTER_inside
!		write(*,*)'rank:',rank, 'COUNTER_overlap', COUNTER_overlap
!	end if 	
!	
!	Tini = Tini + dt
!!	if(rank==0) then
!!		print*, ALL_PARTICLES(1:100)%Xp-Xparticle_read(1,1:100)
!!	end if
!END DO	
	
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
