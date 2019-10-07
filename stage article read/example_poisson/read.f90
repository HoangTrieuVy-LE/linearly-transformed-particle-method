PROGRAM read

  USE mpi

  IMPLICIT NONE

  INTEGER, DIMENSION(MPI_STATUS_SIZE)           :: status
  INTEGER                                       :: rank, code, fh, ntx, nty
  integer, parameter                   :: dp = kind(1.d0)
  REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)    :: u_read
  INTEGER(KIND=MPI_OFFSET_KIND)                  :: file_size
  INTEGER                                        :: double_size
  
  CALL MPI_INIT(code)

  CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, code)

  OPEN(11, FILE='poisson.data', STATUS='OLD')
  READ(11, *) ntx
  READ(11, *) nty
  CLOSE(11)

  ALLOCATE(u_read(ntx, nty))
  u_read(:, :) = 0.d0

  CALL MPI_FILE_OPEN(MPI_COMM_WORLD, "data.dat", &
       MPI_MODE_RDONLY, &
       MPI_INFO_NULL, fh, code)

  !Error checking
  IF (code /= MPI_SUCCESS) THEN
     PRINT *, 'Error opening the file'
     CALL MPI_ABORT(MPI_COMM_WORLD, 2, code)
     CALL MPI_FINALIZE(code)
  END IF

  CALL MPI_FILE_GET_SIZE(fh, file_size, code)
  CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, double_size, code)
  if (file_size /= ntx*nty*double_size) then
    print *, " ATTENTION data.dat does not have the good size (",file_size,",",&
             ntx*nty*double_size,")"
    write(11,*) 0
  else
    CALL MPI_FILE_READ(fh, u_read, SIZE(u_read), &
                       MPI_DOUBLE_PRECISION, status, code)
    WRITE(11, 101)  u_read
101 FORMAT (E19.12)
  end if

  CALL MPI_FILE_CLOSE(fh, code)
  
  CALL MPI_FINALIZE(code)

END PROGRAM read
