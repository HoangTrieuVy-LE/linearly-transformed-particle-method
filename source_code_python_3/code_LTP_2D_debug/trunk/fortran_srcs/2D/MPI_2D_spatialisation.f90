!---------------------------------------------------------------------------  
!> @author 
!> HoangTrieuVy-LE.
!
! DESCRIPTION: 
!> Decomposition of the domain, adapt number of block to number of process, determination of the neighbor.
!--------------------------------------------------------------------------- 

! IMPORT MODULE INSDE SOUBROUTINE, TODO decide which module would be necessary
module MPI_2D_spatialisation_modf90

include 'mpif.h'
use calculsfor_rec_modf90 !include: xpart,mpart
use calculsfor_var_modf90
use calculsfor_ini_modf90



! include 'C:\lib\MPI\Include\mpif.h' ! For windows execution

! Module OPENMP TODO



implicit none 
! DECLARATIONS for MPI
	integer :: rank
	integer :: nb_proc
	integer :: code, i
	integer :: nbdims = 2 !find in config to simply later
	integer, dimension(nbdims): dims
	logical, dimension(nbdims): periods
	integer, dimension(nbdims): coords
	logical :: reorder = .true.


	integer, parameter :: nb_neighbours_2D = 8
	integer, parmeter :: nb_neighbours_3D = 26

! DECLARATION for necessary input data: domain, Xg, Yg, tables of particles, Xpart?, Xcoord?,etc... TODO 

contains
	
	subroutine environnement_initialisation
		call MPI_INIT( code )
		call MPI_COMM_SIZE(MPI_COMM_WORLD,nb_proc,code)
    	call MPI_COMM_RANK(MPI_COMM_WORLD, rank, code)
	end subroutine environnement_initialisation


	subroutine adaptation_block_processor
	! TODO Divide the domain depending on the number of processor
	

	end subroutine_block_processor

	subroutine topology_initialisation*
  		!************
    	!Creation of the Cartesian topology
    	!************
		dims(:) = 0.
		call MPI_DIMS_CREATE(nb_proc, ndims, dims, code)
		periods(:) = .false.
		call MPI_CARTE_CREATE(MPI_COMM_WORLD, ndims, dims, periods, reorder,comm2d,code)
		
	 	if (rank == 0) then
			write (*,'(A)')'------------------------------------------'
			write (*,'(A,i4,A)')'LTP supported by MPI with ',nb_proc,'processes'
			write (*,'(A,i4,A,i4,A)') 'Dimension for the topology: ', &
            						dims(1), ' along x, ', dims(2), ' along  y'			
			write (*,'(A,i4,A,i4)''Domain coords: xmin=, xmax=, ymin, ymax='
			write (*,'(A)')'------------------------------------------'
		end if
	end subroutine topology_initialisation

	subroutine mpi_2d_spatial_decomposition

		!!!!!!!!!!!!!!!!!!!!!!!
		!    iniialisation    !
		!!!!!!!!!!!!!!!!!!!!!!!

		! in 2d, maximum 8 neighbor blocks around 
		rang_left = mpi_proc_null
		rang_right = mpi_proc_null
		rang_up = mpi_proc_null
		rang_down = mpi_proc_null
		rang_up_left = mpi_proc_nulls
		rang_up_right = mpi_proc_null
		rang_down_left = mpi_proc_null
		rang_down_right = mpi_proc_null

		! todo in 3d, maximum 26 neighbor block around
	
		! create the cartesian topology
		! TODO

		call mpi_dims_create(nb_proc,nbdims,dims,code) 
		! dims has ndims size specifying the number of nodes in each dimension


		! todo by now, we do not implement the periodicity of processor to simplify the neighbor find
		! we would like to adapt the number of processor to the number of block, ie we will count the number maximum of processor and assign every block for every processor we have


		periods(:) = .false.

		call mpi_carte_create(mpi_comm_world, nbdims, dims, periods, reorder, comm2d, code)
	
		! calculate the coordinates in the topology
		! todo	
		call mpi_cart_coords(comm2d, rang, nbdims, coords, code)

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!     find neighbors - conditions     !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	
		! neighbor in left_right positions
		if (dim(1)>1) then
			call mpi_cart_shift(comm2d,0,1,rang_left,rang_right,code)
		! todo periodicity?
		end if 
	

		! find neighbors in up_down postion
		if (dim(2)>1) then
			call mpi_cart_shift(comm2d,1,1,rang_down_,rang_up,code)
		! todo periodicity?
	
		end if 

		! find neighbors in the corner in plane x-y
		if (dims(1) > 1 .and. dims(2) > 1) then

		! block has no x=0 and y= dims(2)-1, these blocks can have a up-left neighbor 
			if ( coords(1) .ne. 0 .and. coords(2) .ne. dims(2)-1 ) then
			rang_up_left = rang - (dims(2) - 1) 
			! todo periodicity?

			end if
		! block has no x= dims(1)-1 and y = dims(2)-1, these blocks can have a up-right neighbor
			if ( coords(1) .ne. dim(1)-1 .and. coords(2) .ne. dims(2)-1) then
			rang_up_right = rang + dims(2) + 1 
			! todo periodicity?
			end if

		! block has no x = 0 and y = 0, these blocks can have a down-left neighbor
			if (coords(1) .ne. 0 .and. coords(2) .ne. 0) then 
				rang_down_left = rang - dims(2) -1
			! todo periodicity?		
			end if

		! block has no x = dims(1)-1 and y =0, these blocks can have a down-right neighbor
			if (coords(1) .ne. dims(1)-1 and coords(2) .ne. 0) then
				rand_down_right = rang + dims(2) -1
			! todo periodicity?	
			end if

		end if ! end if for neighbor in plane x-y


		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!  todo add general case dims(3)>1  !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	end subroutine mpi_2d_spatial_decomposition

	subroutine 	environnement_finalization
		call MPI_FINALIZE(code)
	end subroutine environnement_finalization

end module MPI_2D_spatialisation_modf90
