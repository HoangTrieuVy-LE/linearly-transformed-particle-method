
PROGRAM hello_world_mpi
INCLUDE 'C:\lib\MPI\Include\mpif.h'

integer rank,size,ierror

call MPI_INIT(ierror)
call MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierror)
call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)

! DO i=0,3,1
! 	IF(i==rank) THEN
! 	print*, "From process", rank, "of ", size
! 	END IF
! 	print*,i
! 	call MPI_BARRIER(MPI_COMM_WORLD,ierror)
! END DO

integer :: i,rank,ntasks,compte,init,fin,nloops

call MPI_INIT(i, rank, ntasks, compte,init,fin, nloops, total_nloops)

call MPI_COMM_RANK(MPI_COMM_WORLD,rank)
call MPI_COMM_SIZE(MPI_COMM_WORLD,ntasks)



compte = 250
init =rank*compte
fin = init+compte

nloops = 0
DO i =  init, fin-1,1
	nloops = nloops + 1
END DO

print*, "Task",rank,"performed",nloops


call MPI_FINALIZE(ierror)

END PROGRAM