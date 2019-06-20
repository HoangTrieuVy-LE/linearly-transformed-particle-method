!module_parallel_mpi.f90
!!!!
!!subroutine env_init
!!subroutine topology_init
!!subroutine domain_boundaries
!!subroutine domain_neighbours
!!subroutine derived_datatypes
!!subroutine communications
!!subroutine nonblocking_communications
!!function   global_error
!!subroutine mpi_write
!!subroutine env_finalize
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE parallel
  USE TYPE_PARAMS
  USE MPI
  IMPLICIT NONE

  !Local Sub-Domain rank
  INTEGER                                   :: rank
  !Number of processes
  INTEGER                                   :: size
  !Communicator of the Cartesian topology
  INTEGER                                   :: comm2d
  !Number of dimensions for the grid
  INTEGER, PARAMETER                        :: ndims = 2
  !Number of processes in each dimension for the Cartesian topology
  INTEGER, DIMENSION(ndims)                 :: dims
  !Topology periodicity
  LOGICAL, DIMENSION(ndims)                 :: periods
  !Coordinates of the local domain
  INTEGER, DIMENSION(ndims)                 :: coords
  !Array storing the rank of neighbours
  INTEGER, PARAMETER                        :: NB_NEIGHBOURS = 4
  INTEGER, PARAMETER                        :: N=1, E=2, S=3, W=4
  INTEGER, DIMENSION(NB_NEIGHBOURS)         :: neighbour
  !Derived datatypes
  INTEGER                                   :: type_line, type_column
  !MPI
  INTEGER                                   :: code  

CONTAINS

  SUBROUTINE env_init
    !************
    !Initialization of the MPI environnement
    !************

    !MPI initialization
    CALL MPI_INIT(code)

    !Who I am
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, code)

    !Total number of processes
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, size, code)

  END SUBROUTINE env_init

  SUBROUTINE topology_init
    !************
    !Creation of the Cartesian topology
    !************

    !MPI
    LOGICAL, PARAMETER                        :: reorganisation = .FALSE.

    !Read ntx and nty in the file poisson.data
    OPEN(10, FILE='poisson.data', STATUS='OLD')
    READ(10, *) ntx
    READ(10, *) nty
    CLOSE(10)

    !Number of processes on each dimension (depends on the total number of processes)
    dims(:)= 0
    CALL MPI_DIMS_CREATE(size, ndims, dims, code)

    !Creation of the 2D cartesian topology (no periodicity)
    periods(:) = .FALSE.
    CALL MPI_CART_CREATE(MPI_COMM_WORLD, ndims, dims, periods, &
         reorganisation, comm2d, code)

    IF (rank == 0) THEN
       WRITE (*,'(A)') '-----------------------------------------'
       WRITE (*,'(A,i4,A)') 'Execution poisson with ', size, ' MPI processes'
       WRITE (*,'(A,i4,A,i4)') 'Size of the domain: ntx=', ntx, ' nty=', nty
       WRITE (*,'(A,i4,A,i4,A)') 'Dimension for the topology: ', &
            dims(1), ' along x, ', dims(2), ' along  y'
       WRITE (*,'(A)') '-----------------------------------------'
    END IF

  END SUBROUTINE topology_init


  SUBROUTINE domain_boundaries
    !************
    !Computation of the local grid boundary coordinates (global indexes)
    !************

    ! Coordinates
    CALL MPI_CART_COORDS(comm2d, rank, ndims, coords, code)

    ! X-axis limits
    sx = (coords(1)*ntx)/dims(1)+1
    ex = ((coords(1)+1)*ntx)/dims(1)

    ! Y-axis limits
    sy = (coords(2)*nty)/dims(2)+1
    ey = ((coords(2)+1)*nty)/dims(2)

    WRITE (*,'(A,i4,A,i4,A,i4,A,i4,A,i4,A)') 'Rank in the topology: ', rank, &
         ' Local Grid Index:',  sx, ' to', ex, ' along x, ', &
         sy, ' to', ey, ' along y'

  END SUBROUTINE domain_boundaries

  SUBROUTINE domain_neighbours
    !************
    !Neighbours 
    !************

    !Initialization of the neighbor array
    neighbour(:) = MPI_PROC_NULL

    !Get my northern and southern neighbours
    CALL MPI_CART_SHIFT(comm2d, 0, 1, neighbour(N), neighbour(S), code)

    !Get my western and eastern neighbours
    CALL MPI_CART_SHIFT(comm2d, 1, 1, neighbour(W), neighbour(E), code)

    WRITE (*,'(A,i4,A,i4,A,i4,A,i4,A,i4)') "Process ", rank, " neighbour: N", neighbour(N), " E", neighbour(E), &
            " S", neighbour(S), " W", neighbour(W)

  END SUBROUTINE domain_neighbours


  SUBROUTINE derived_datatypes
    !************
    !Creation of the derived datatypes needed to exchange points with neighbours
    !************

    !Creation of the type_line derived datatype to exchange points
    !with northern to southern neighbours
    CALL MPI_TYPE_VECTOR(ey-sy+1, 1, ex-sx+3, &
         MPI_DOUBLE_PRECISION, type_line, code)
    CALL MPI_TYPE_COMMIT(type_line, code)

    !Creation of the type_column derived datatype to exchange points
    !with western to eastern neighbours
    CALL MPI_TYPE_CONTIGUOUS(ex - sx + 1, MPI_DOUBLE_PRECISION, &
         type_column, code)
    CALL MPI_TYPE_COMMIT(type_column, code)

  END SUBROUTINE derived_datatypes


  SUBROUTINE communications(u)
    !************
    !Exchange the points at the interface
    !************

    REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :), INTENT(inout) :: u

    !MPI constants
    INTEGER, PARAMETER                   :: tag=100
    INTEGER, DIMENSION(MPI_STATUS_SIZE)  :: status

    !Send to neighbour N and receive from neighbour S
    CALL MPI_SENDRECV(u(sx, sy), 1,   type_line,           neighbour(N), &
         tag,  u(ex+1, sy), 1,        type_line,           neighbour(S), &
         tag, comm2d, status, code)

    !Send to neighbour S and receive from neighbour N
    CALL MPI_SENDRECV(u(ex, sy), 1,   type_line,           neighbour(S), &
         tag,  u(sx-1, sy), 1,        type_line,           neighbour(N), &
         tag, comm2d, status, code)

    !Send to neighbour W  and receive from neighbour E
    CALL MPI_SENDRECV(u(sx, sy), 1, type_column,           neighbour(W), &
         tag,  u(sx, ey+1), 1, type_column,                neighbour(E), &
         tag, comm2d, status, code)

    !Send to neighbour E  and receive from neighbour W
    CALL MPI_SENDRECV(u(sx, ey), 1, type_column,           neighbour(E), &
         tag,  u(sx, sy-1), 1, type_column,                neighbour(W), &
         tag, comm2d, status, code)

  END SUBROUTINE communications

  FUNCTION global_error(u, u_new)
    !************
    !Calcul for the global error (maximum of the locals errors)
    !************
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :), INTENT(in) :: u, u_new

    REAL(kind=dp)              :: global_error, local_error

    local_error = MAXVAL (ABS(u(sx:ex, sy:ey) &
                             -u_new(sx:ex, sy:ey)))

    !Calcul of local error on each sub-domain
    CALL MPI_AllREDUCE(local_error, global_error, 1, MPI_DOUBLE_PRECISION, &
         MPI_MAX, comm2d, code)
  END FUNCTION global_error

  SUBROUTINE mpi_write(u)
    !********************
    ! Write array u inside a domain for each process in the data.dat file
    !********************
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :), INTENT(inout) :: u

    INTEGER, DIMENSION(MPI_STATUS_SIZE)   :: status
    INTEGER                               :: fh
    INTEGER(kind = MPI_OFFSET_KIND)       :: initial_displacement
    INTEGER, PARAMETER                    :: array_rank=2
    INTEGER, DIMENSION(array_rank)        :: shape_array, shape_sub_array, start_coord
    INTEGER, DIMENSION(array_rank)        :: shape_view_array, shape_sub_view_array, start_view_coord
    INTEGER                               :: type_sub_array, type_sub_view_array

    !
    !Open file "data.dat" in write mode
    !
    CALL MPI_FILE_OPEN(comm2d, "data.dat", &
         MPI_MODE_WRONLY + MPI_MODE_CREATE, &
         MPI_INFO_NULL, fh, code)

    !Error checking
    IF (code /= MPI_SUCCESS) THEN
       PRINT *, 'Error opening the file'
       CALL MPI_ABORT(comm2d, 2, code)
    END IF

    !
    ! Creation of the derived datatype type_sub_array corresponding to the matrix u without ghost cells
    !
    !Shape  of the array 
    shape_array(:)= SHAPE(u)

    !Shape of the subarray
    shape_sub_array(:) = SHAPE(u(sx:ex, sy:ey))

    !Starting coordinates of the subarray
    start_coord(:) = (/ 1 , 1  /)

    !Creation of the derived datatype type_sub_array
    CALL MPI_TYPE_CREATE_SUBARRAY(array_rank, shape_array, shape_sub_array, start_coord, &
         MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, type_sub_array, code)

    !Commit type_sub_array
    CALL MPI_TYPE_COMMIT(type_sub_array, code)

    !
    !Creation of the derived datatype type_sub_view_array for the view on the file
    !
    !Shape of the array
    shape_view_array(:)= (/ ntx, nty /)

    !Shape of the subarray
    shape_sub_view_array(:) = SHAPE(u(sx:ex, sy:ey))

    !Starting coordinates of the subarray
    start_view_coord(:) =  (/ sx-1 , sy-1 /)

    !Creation of the derived datatype type_sub_view_array
    CALL MPI_TYPE_CREATE_SUBARRAY(array_rank, shape_view_array, shape_sub_view_array, start_view_coord, &
         MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, type_sub_view_array, code)

    !Commit type_sub_view_array
    CALL MPI_TYPE_COMMIT(type_sub_view_array, code)

    !
    !Change the file view
    !
    initial_displacement = 0
    CALL MPI_FILE_SET_VIEW(fh, initial_displacement, MPI_DOUBLE_PRECISION, &
         type_sub_view_array, "native", MPI_INFO_NULL, code)

    !
    !Write u for each process with the view
    !
    CALL MPI_FILE_WRITE_ALL(fh, u, 1, type_sub_array, status, code)

    !
    !Close file
    !
    CALL MPI_FILE_CLOSE(fh, code)

  END SUBROUTINE mpi_write

  SUBROUTINE env_finalize
    !************
    !Terminates MPI execution environment
    !************

    CALL MPI_FINALIZE(code)

  END SUBROUTINE env_finalize

END MODULE parallel
