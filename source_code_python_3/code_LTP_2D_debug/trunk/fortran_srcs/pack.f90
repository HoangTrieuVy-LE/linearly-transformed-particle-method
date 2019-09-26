!---------------------------------------------------------------------------  
!> @author 
!> HoangTrieuVy-LE.
!
! DESCRIPTION: 
!> Name: PACKMPI
!> Stored global variables, the detail informations are below
!--------------------------------------------------------------------------- 

MODULE PACKMPI

USE MPI

IMPLICIT NONE

	! DECLARATIONS for MPI
	
	! Error variable for MPI 
	INTEGER                      :: code

	! Core rank
	INTEGER                      :: rank
	
	! Number of proc
	INTEGER                      :: nb_proc

	! Loop index
	INTEGER                      :: i,j


	! COMMUNICATOR OF THE CARTESIAN TOPOLOGY
  	INTEGER                      :: comm2d

	! number of dimension
	INTEGER, PARAMETER           :: nbdims = 2 

	! Indicate the number of cores for each dimension
	INTEGER, DIMENSION(nbdims)   :: dims
	
	! In case there is a period of core but we do not make here
	LOGICAL, DIMENSION(nbdims)   :: periods

	! Topology Coordinates in MPI of each proc
	INTEGER, DIMENSION(nbdims)   :: coords

	LOGICAL                      :: reorder = .true.

	! The lenght size of each block, and the coordinates in simulation region for each block  
	DOUBLE PRECISION             :: block_step_x, block_step_y, start_x, start_y, end_x, end_y
	
	! Neighbour rank
	INTEGER                      :: rank_left, rank_right, rank_up, rank_down,&
 									rank_up_left, rank_up_right, rank_down_left,&
								 	rank_down_right
	! Number of neighbour proc in 2D
	INTEGER, PARAMETER           :: NB_NEIGHBOURS = 8
	INTEGER			             :: nb_neighbours_actually
	
	!INTEGER(kind=MPI_OFFSET_KIND)                 :: offset


END MODULE PACKMPI

