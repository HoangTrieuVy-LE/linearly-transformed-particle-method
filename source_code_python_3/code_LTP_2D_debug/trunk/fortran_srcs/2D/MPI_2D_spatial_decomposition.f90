!---------------------------------------------------------------------------  
!> @author 
!> HoangTrieuVy-LE.
!
! DESCRIPTION: 
!> Decomposition of the domain, adapt number of block to number of process, determination of the neighbor.
!--------------------------------------------------------------------------- 

! IMPORT MODULE INSDE SOUBROUTINE, TODO decide which module is necessary
use calculsfor_rec_modf90
use calculsfor_var_modf90
use calculsfor_ini_modf90



subroutine MPI_2D_spatial_decomposition
	! DECLARATIONS
	! TODO

	! INIIALISATION
	! TODO
	! En 2D, maximum 8 neighbor blocks around 
	rang_left = MPI_PROC_NULL
	rang_right = MPI_PROC_NULL
	rang_up = MPI_PROC_NULL
	rang_down = MPI_PROC_NULL
	rang_up_left = MPI_PROC_NULL
	rang_up_right = MPI_PROC_NULL
	rang_down_left = MPI_PROC_NULL
	rang_down_right = MPI_PROC_NULL
	
	! Create the cartesian topology
	! TODO
	! Calculate the coordinantes in the topology
	! TODO	
	! Find neighbors - conditions 
	! TODO


end subroutine MPI_2D_spatial_decomposition


