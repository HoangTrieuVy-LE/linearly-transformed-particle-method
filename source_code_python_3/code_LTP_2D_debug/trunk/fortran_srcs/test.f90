program test

include 'mpif.h'
integer :: size, code
call MPI_INIT(code)


call MPI_COMM_SIZE(MPI_COMM_WORLD,size, code)
print*, size
call MPI_FINALIZE(code)
end program test
