program test
implicit none

include 'mpif.h'

integer,parameter:: nb_lines=5,nb_columns=6
integer,parameter:: tag=100
real,dimension(nb_lines,nb_columns)   :: a
integer,dimension(MPI_STATUS_SIZE)    :: msgstatus
integer:: rank,code,type_column

call MPI_INIT(code)
call MPI_COMM_RANK(MPI_COMM_WORLD,rank,code)
a(:,:) =real(rank)

call MPI_TYPE_CONTIGUOUS(nb_lines,MPI_REAL,type_column,code)
call MPI_TYPE_COMMIT(type_column,code)

if( rank == 0 )then
call MPI_SEND(a(1,1),1,type_column,1,tag,MPI_COMM_WORLD,code)
elseif( rank == 1 )then
call MPI_RECV(a(1,nb_columns),nb_lines,MPI_REAL,0,tag,&
	MPI_COMM_WORLD,msgstatus,code)
end if
call MPI_TYPE_FREE(type_column,code)
call MPI_FINALIZE(code)
call MPI_FINALIZE(code)

end program test
