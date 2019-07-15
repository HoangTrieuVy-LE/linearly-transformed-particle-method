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

INTEGER   :: Nboucle, Npic


Nboucle = 0
Npic    = 0


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

!	!-------------------------------------------------------------------!
!	!!!  PARTICLES ATTRIBUTION  &  OVERLAP PARTICLES LISTS CREATION   !!!
!	!-------------------------------------------------------------------!

CALL initiation_table
CALL parttype_convert
CALL neighbouring

CALL particle_distribution


!if(rank==4) then
!	print*,'FJ=1,BJ=2,FK=3,BK=4,FJBK=5,BJFK=6,BJBK=7,FJFK=8'
!	print*,'rank',rank,'  COUNTER-INSIDE',COUNTER_inside
!	DO i=1,8
!!	print*,'Danger on',i,':',COUNTER_danger(i,rank)
!	print*,'Overlap on',i,':',COUNTER_overlap(i,rank)
!	END DO
!end if

!	!-------------------------------------------------------------------!
!	!!!LOOP ON PARTICLES INSIDE BLOCK AND THE OTHERS IN OVERLAP TABLES!!!
!	!-------------------------------------------------------------------!
 if(rank==4) then
print*,'----------------------------------'
	write(*,*)'rank:',rank, 'COUNTER_inside', COUNTER_inside(rank)
	write(*,*)'              				   UP	','	DOWN	','  RIGHT','       LEFT	',' UP-LEFT  ',' DOWN-RIGHT','    DOWN-LEFT','  UP-RIGHT'	
	write(*,*)'rank:',rank, 'COUNTER_overlap', COUNTER_overlap(:,rank)
	write(*,*)'rank:',rank, 'COUNTER_danger', COUNTER_danger(:,rank)
end if
DO WHILE(T_start<=T_end)


	CALL send_overlap_and_danger_particle

	CALL block_loop_on_block_global_table

	!-------------------------------------------------------------------!
	!!!              UPDATE ALL PARTICLES INFORMATIONS                !!!
	!-------------------------------------------------------------------!

	CALL update_ALL_particles
	


	T_start = T_start + time_step

	Npic = Npic + 1
	Nboucle = Nboucle + 1

END DO	
 if(rank==4) then
print*,'----------------------------------'
	write(*,*)'rank:',rank, 'COUNTER_inside', COUNTER_inside(rank)
	write(*,*)'              				   UP	','	DOWN	','  RIGHT','       LEFT	',' UP-LEFT  ',' DOWN-RIGHT','    DOWN-LEFT','  UP-RIGHT'	
	write(*,*)'rank:',rank, 'COUNTER_overlap', COUNTER_overlap(:,rank)
	write(*,*)'rank:',rank, 'COUNTER_danger', COUNTER_danger(:,rank)
end if

	!-------------------------------------------------------------------!
	!!!                        UPDATE FILE IN                         !!!
	!-------------------------------------------------------------------!
CALL update_all_particle_information


OPEN(unit= 4, FILE='trunk/fortran_srcs/temp_out.txt')			
	write(4,5) T_start-time_step
	write(4,6) Npic
	write(4,7) Nboucle
	5 	format (f16.10)
	6	format (10i7)
	7	format (10i7)
	close(4)



	!-------------------------------------------------------------------!
	!!!                         DEALLOCATION                          !!!
	!-------------------------------------------------------------------!
		
!	CALL dealloc_XMD
!	CALL dealloc_all_particle_table


CALL environnement_finalization
END PROGRAM tri_casier_modf90
