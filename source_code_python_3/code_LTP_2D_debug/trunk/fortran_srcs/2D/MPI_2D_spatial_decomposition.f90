!---------------------------------------------------------------------------  
!> @author 
!> HoangTrieuVy-LE.
!
! DESCRIPTION: 
!> Decomposition of the domain, adapt number of block to number of process, determination of the neighbor.
!--------------------------------------------------------------------------- 
subroutine MPI_2D_spatial_decomposition
	! DECLARATIONS
	

	! INIIALISATION
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
	! Calculate the coordinantes in the topology
	! Find neighbors - conditions TODO


end subroutine MPI_2D_spatial_decomposition


