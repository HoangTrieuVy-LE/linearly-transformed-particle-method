!---------------------------------------------------------------------------  
!> @author 
!> HoangTrieuVy-LE.
!
! DESCRIPTION: 
!> Decomposition of the domain, adapt number of block to number of process, determination of the neighbour.
!--------------------------------------------------------------------------- 

! IMPORT MODULE INSDE SOUBROUTINE, TODO decide which module would be necessary
MODULE mpi_2D_spatialisation_modf90

IMPLICIT NONE 


INCLUDE 'mpif.h'
USE mod_particle_2D_modf90

	! DECLARATIONS for MPI
	INTEGER                      :: rank

	INTEGER                      :: nb_proc

	INTEGER                      :: code, i

	! COMMUNICATOR OF THE CARTESIAN TOPOLOGY
  	INTEGER                      :: comm2d

	INTEGER, PARAMETER           :: nbdims = 2 !find in config to simply later

	INTEGER, DIMENSION(nbdims)   :: dims

	LOGICAL, DIMENSION(nbdims)   :: periods

	INTEGER, DIMENSION(nbdims)   :: coords

	LOGICAL                      :: reorder = .true.

	DOUBLE PRECISION             :: block_step_x, block_step_y, start_x, start_y, end_x, end_y

	INTEGER, PARAMETER           :: nb_neighbours_2D = 8

	INTEGER, PARAMETER           :: nb_neighbours_3D = 26

	INTEGER                      :: rank_left, rank_right, rank_up, rank_down,&
 									rank_up_left, rank_up_right, rank_down_left,&
								 	rank_down_right



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
		CALL MPI_DIMS_CREATE(nb_proc, nbdims, dims, code)

		! Creation of the 2D cartesian topology (no periodicity)
		periods(:) = .false.
		CALL MPI_CART_CREATE(MPI_COMM_WORLD, nbdims, dims, periods, reorder,comm2d,code)
		
		IF (rank == 0) THEN
			WRITE (*,'(A,i4,A,i4,A)') 'Dimension for the topology: ', &
							dims(1), ' along x, ', dims(2), ' along  y'
		END IF

	END SUBROUTINE topology_initialisation


	SUBROUTINE set_block_grid(Lx1,Lx2,Ly1,Ly2)
		DOUBLE PRECISION, INTENT(in) :: Lx1,Lx2,Ly1,Ly2		
		
		! TODO Divide the domain depending on the number of processor
		! Consider we have Lx1,Lx2,Ly1,Ly2
		! After having dims(), we can find block_step_x and block_step_y
		block_step_x = (Lx2-Lx1)/dims(1)
		block_step_y = (Ly2-Ly1)/dims(2)

		! calculate the coordinates in the topology
		CALL mpi_cart_coords(comm2d, rank, nbdims, coords, code)		
		
		! Limits at X-axis
		start_x = (coords(1)*block_step_x)/dims(1)
		end_x =((coords(1)+1)*block_step_x)/dims(1)

		! Limits at Y-axis
		start_y = (coords(2)*block_step_y)/dims(2)
		end_y = ((coords(2)+1)*block_step_y)/dims(2)

		WRITE (*,'(A,i4,A,F5.2,A,F5.2,A,F5.2,A,F5.2,A)') 'Rank in the topology: ', rank, &
         ' Local Grid Index:',  start_x, ' to', end_x, ' along x, ', &
         start_y, ' to', end_y, ' along y'

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
		!     find neighbours - conditions     !
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

	SUBROUTINE 	environnement_finalization
		CALL MPI_FINALIZE(code)
	END SUBROUTINE environnement_finalization

END MODULE mpi_2D_spatialisation_modf90
