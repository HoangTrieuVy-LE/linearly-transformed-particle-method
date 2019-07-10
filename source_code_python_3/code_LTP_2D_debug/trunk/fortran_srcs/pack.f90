MODULE PACKMPI

USE MPI

IMPLICIT NONE

	! DECLARATIONS for MPI
	INTEGER                      :: code

	INTEGER                      :: rank
	
	INTEGER                      :: nb_proc

	INTEGER                      :: i

	! COMMUNICATOR OF THE CARTESIAN TOPOLOGY
  	INTEGER                      :: comm2d

	INTEGER, PARAMETER           :: nbdims = 2 

	INTEGER, DIMENSION(nbdims)   :: dims

	LOGICAL, DIMENSION(nbdims)   :: periods

	INTEGER, DIMENSION(nbdims)   :: coords

	LOGICAL                      :: reorder = .true.

	DOUBLE PRECISION             :: block_step_x, block_step_y, start_x, start_y, end_x, end_y

	INTEGER                      :: rank_left, rank_right, rank_up, rank_down,&
 									rank_up_left, rank_up_right, rank_down_left,&
								 	rank_down_right

	INTEGER, PARAMETER           :: NB_NEIGHBOURS = 8
	INTEGER			             :: nb_neighbours_actually
	
	!INTEGER(kind=MPI_OFFSET_KIND)                 :: offset


END MODULE PACKMPI

