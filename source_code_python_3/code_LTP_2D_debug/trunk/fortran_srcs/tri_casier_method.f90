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

Npic =0
Nboucle = 0

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
!	if (rank==4) then
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

!print*,'rank',rank,coords(1),coords(2),'has:',start_x,end_x,'and',start_y,end_y
!	!-------------------------------------------------------------------!
!	!!!  PARTICLES ATTRIBUTION  &  OVERLAP PARTICLES LISTS CREATION   !!!
!	!-------------------------------------------------------------------!

CALL initiation_table
CALL parttype_convert
CALL neighbouring


CALL neighbour_limit_finding

CALL particle_distribution_v2

!	!-------------------------------------------------------------------!
!	!!!LOOP ON PARTICLES INSIDE BLOCK AND THE OTHERS IN OVERLAP TABLES!!!
!	!-------------------------------------------------------------------!

if(rank==6.or.rank==7.or.rank==8.or.rank==11.or.rank==12.or.rank==13.or.rank==16.or.rank==17.or.rank==18) then
print*,'----------------------------------'
	write(*,*)'rank:',rank, 'COUNTER_inside', COUNTER_inside(rank)
	write(*,*)'              				   UP	','	DOWN	','  RIGHT','       LEFT	',' UP-LEFT  ',' DOWN-RIGHT','    DOWN-LEFT','  UP-RIGHT'	
	write(*,*)'rank:',rank, 'COUNTER_overlap', COUNTER_overlap(:,rank)
	write(*,*)'rank:',rank, 'COUNTER_danger', COUNTER_danger(:,rank)
!print*,'IND inside',IND_inside(1:4)

end if


if(rank==8) then
!DO i=1,number_of_particles
!		print*,ALL_PARTICLES(i)%ID,ALL_PARTICLES(i)%Xp,ALL_PARTICLES(i)%Yp
!	END DO
!DO i=1,8
!	write(*,*)'rank:',rank,neighbour_limit(:,i)
!END DO

end if


CALL send_overlap_and_danger_particle

!CALL block_loop_on_block_global_table


!CALL update_ALL_particles

!CALL update_all_particle_information


DO WHILE(T_start<=T_end)

	CALL block_loop_on_block_global_table
!	
!	!-------------------------------------------------------------------!
!	!!!              UPDATE ALL PARTICLES INFORMATIONS                !!!
!	!-------------------------------------------------------------------!

	CALL update_ALL_particles

!if(rank==5.or.rank==6.or.rank==9.or.rank==10) then
!print*,'STEP: ',Npic
!print*,'----------------------------------'
!	write(*,*)'rank:',rank, 'COUNTER_inside', COUNTER_inside(rank)
!	write(*,*)'              				   UP	','	DOWN	','  RIGHT','       LEFT	',' UP-LEFT  ',' DOWN-RIGHT','    DOWN-LEFT','  UP-RIGHT'	
!	write(*,*)'rank:',rank, 'COUNTER_overlap', COUNTER_overlap(:,rank)
!	write(*,*)'rank:',rank, 'COUNTER_danger', COUNTER_danger(:,rank)
!write(*,*)'rank:',rank, 'COUNTER_leave', COUNTER_leave(:,rank)
!end if



	T_start = T_start + time_step

	Npic = Npic + 1
	Nboucle = Nboucle + 1

END DO	

	!-------------------------------------------------------------------!
	!!!                        UPDATE FILE IN                         !!!
	!-------------------------------------------------------------------!
CALL update_all_particle_information


OPEN(unit= 4, FILE='trunk/fortran_srcs/temp_out.txt')			
	write(4,5) T_start
	write(4,6) Npic
	write(4,7) Nboucle
	5 	format (f16.10)
	6	format (10i7)
	7	format (10i7)
	close(4)



!	!-------------------------------------------------------------------!
!	!!!                         DEALLOCATION                          !!!
!	!-------------------------------------------------------------------!
!		
!	CALL dealloc_XMD
	CALL dealloc_all_particle_table


CALL environnement_finalization
END PROGRAM tri_casier_modf90
