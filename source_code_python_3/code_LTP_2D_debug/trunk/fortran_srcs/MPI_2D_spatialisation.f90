!---------------------------------------------------------------------------  
!> @author 
!> HoangTrieuVy-LE.
!
! DESCRIPTION: 
!> Decomposition of the domain, adapt number of block to number of process, determination of the neighbour.
!--------------------------------------------------------------------------- 

! IMPORT MODULE INSDE SOUBROUTINE, TODO decide which module would be necessary
MODULE mpi_2D_spatialisation_modf90

USE data_launch
USE mod_particle_2D_modf90

IMPLICIT NONE 


CONTAINS
	
	SUBROUTINE environnement_initialisation
		CALL MPI_INIT( code )

		CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nb_proc,code)

    	CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, code)
	END SUBROUTINE environnement_initialisation

		! By now, we do not implement the periodicity of processor to simplify the neighbour find
		! We would like to adapt the number of processor to the number of block, ie we will count the number of processor and assign every block for every processor we have

	SUBROUTINE topology_initialisation

  		!**********************************!
    	!Creation of the Cartesian topology!
    	!**********************************!

		! Number of processes on each dimension (depends on the total number of processes)
		dims(:) = 0

!       Example for dims(1) = 3 and dims(2) = 3
!		-------------
!		! 2 ! 5 ! 8 !
!		-------------
!		! 1 ! 4 ! 7 !
!		-------------
!		! 0 ! 3 ! 6 !
! 		-------------
		

		CALL MPI_DIMS_CREATE(nb_proc, nbdims, dims, code)
		! Creation of the 2D cartesian topology (no periodicity)
		periods(:) = .false.
		CALL MPI_CART_CREATE(MPI_COMM_WORLD, nbdims, dims, periods, reorder,comm2d,code)

	END SUBROUTINE topology_initialisation


	SUBROUTINE set_block_grid(Lx1,Lx2,Ly1,Ly2)
		DOUBLE PRECISION, INTENT(in) :: Lx1,Lx2,Ly1,Ly2		
		INTEGER   :: STATUS=0		
		

!		print*,dims(1)
!		print*,Lx2,Lx1
!		print*,'block-step',block_step_x
!		print*, rank, start_x,end_x

		! calculate the block coordinates in the topology
		CALL mpi_cart_coords(comm2d, rank, nbdims, coords, code)

!                        !       !       !       !         
!		                 !       !       !       !         
!                 (0,4)	 ! (1,4) ! (2,4) ! (3,4) !   (4,4)      
!		                 !       !       !       !         
!                -------------------------------------------
!                 (0,3)	 ! (1,3) ! (2,3) ! (3,3) !  (4,3)      
!                -------------------------------------------
!                 (0,2)  ! (1,2) ! (2,2) ! (3,2) !  (4,2)       
!                -------------------------------------------
!                 (0,1)  ! (1,1) ! (2,1) ! (3,1) !  (4,1)       
!                -------------------------------------------
!	                 	 !       !       !       !         
!                 (0.0)	 ! (1,0) ! (2,0) ! (3,0) !  (4,0)        
!	                 	 !		 !	     !       !         
!                        !       !       !       !         


! LES BLOCK SURROUNDING HAVE INFINITY LENGHT

		if(dims(1) <=2 .or. dims(2)<=2) then
			if(rank==0) then
			print*, "ERROR DIMENSION, MUSTE BE AS LEAST 3x3"
			end if
			call EXIT(STATUS)
		else 
		! After having dims(), we can find block_step_x and block_step_y
		block_step_x = (Lx2-Lx1)/(dims(1)-2)
		block_step_y = (Ly2-Ly1)/(dims(2)-2)

!		-------------------------
!		! (0,2) ! (1,2) ! (2,2) !
!		-------------------------
!		! (0,1) ! (1,1) ! (2,1) !
!		-------------------------
!		! (0,0) ! (1,0) ! (2,0) !
! 		-------------------------
		
		! Limits at X-axis
		start_x = ((coords(1)-1)*block_step_x)+Lx1
		end_x =(coords(1)*block_step_x)+Lx1

		! Limits at Y-axis
		start_y = ((coords(2)-1)*block_step_y)+Ly1
		end_y = (coords(2)*block_step_y)+Ly1

		
		if(coords(1)==0) then
			start_x = -1000000
		end if
		if(coords(1)==dims(1)-1)then
			end_x = 1000000
		end if
		if(coords(2)==0) then
			start_y = -1000000
		end if	
		if(coords(2)==dims(2)-1) then
			end_y = 1000000
		end if
		 
		end if

!		print*,'rank',coords(1),coords(2),'has:',start_x,end_x,'and',start_y,end_y


	END SUBROUTINE set_block_grid

	SUBROUTINE neighbour_blocks
		!!!!!!!!!!!!!!!!!!!!!!!
		!    iniialisation    !
		!!!!!!!!!!!!!!!!!!!!!!!

		! in 2d, maximum 8 neighbour blocks around 
		rank_left = mpi_proc_null
		rank_right = mpi_proc_null
		rank_up = mpi_proc_null
		rank_down = mpi_proc_null
		rank_up_left = mpi_proc_null
		rank_up_right = mpi_proc_null
		rank_down_left = mpi_proc_null
		rank_down_right = mpi_proc_null

		! TODO in 3d, maximum 26 neighbour block around
	
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!     find neighbours - conditions    !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
		! neighbour in left_right positions
		IF (dims(1)>1) THEN
			CALL mpi_cart_shift(comm2d,0,1,rank_left,rank_right,code)
		! todo periodicity?
		END IF 
	

		! find neighbours in up_down postion
		IF (dims(2)>1) THEN
			CALL mpi_cart_shift(comm2d,1,1,rank_down,rank_up,code)
		! todo periodicity?
	
		END IF 

		! find neighbours in the corner in plane 2D
		IF (dims(1) > 1 .and. dims(2) > 1) THEN

		! block has no x=0 and y= dims(2)-1, these blocks can have a up-left neighbour 
			IF ( coords(1) .ne. 0 .and. coords(2) .ne. dims(2)-1 ) THEN
			rank_up_left = rank - (dims(2) - 1) 
			! todo periodicity?

			END IF
		! block has no x= dims(1)-1 and y = dims(2)-1, these blocks can have a up-right neighbour
			IF ( coords(1) .ne. dims(1)-1 .and. coords(2) .ne. dims(2)-1) THEN
			rank_up_right = rank + dims(2) + 1 
			! todo periodicity?
			END IF

		! block has no x = 0 and y = 0, these blocks can have a down-left neighbour
			IF (coords(1) .ne. 0 .and. coords(2) .ne. 0) THEN 
				rank_down_left = rank - dims(2) -1
			! todo periodicity?		
			END IF

		! block has no x = dims(1)-1 and y =0, these blocks can have a down-right neighbour
			IF (coords(1) .ne. dims(1)-1 .and. coords(2) .ne. 0) THEN
				rank_down_right = rank + dims(2) -1
			! todo periodicity?	
			END IF

		END IF ! end if for neighbour in plane x-y


		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!  TODO add general case dims(3)>1  !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	END SUBROUTINE neighbour_blocks


	SUBROUTINE neighbour_counter
		nb_neighbours_actually = 0
		if (rank_left .ge. 0) then 
			nb_neighbours_actually = nb_neighbours_actually+ 1
		end if
		if (rank_right .ge. 0) then 
		nb_neighbours_actually =nb_neighbours_actually+ 1
		end if
		if (rank_up .ge. 0)then
			nb_neighbours_actually =nb_neighbours_actually + 1
		end if
		if (rank_down .ge. 0)then
			nb_neighbours_actually = nb_neighbours_actually+ 1
		end if
		if (rank_up_left .ge. 0)then
			nb_neighbours_actually = nb_neighbours_actually+1
		end if
		if (rank_up_right .ge. 0)then
			nb_neighbours_actually =nb_neighbours_actually+ 1
		end if
		if (rank_down_left .ge. 0)then 
			nb_neighbours_actually =nb_neighbours_actually+ 1
		end if
		if (rank_down_right  .ge. 0) then 
			nb_neighbours_actually =nb_neighbours_actually+ 1
		end if		
	END SUBROUTINE neighbour_counter



	SUBROUTINE 	environnement_finalization
		CALL MPI_FINALIZE(code)
	END SUBROUTINE environnement_finalization

END MODULE mpi_2D_spatialisation_modf90
