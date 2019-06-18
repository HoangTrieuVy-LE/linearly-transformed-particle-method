program TP1 
	implicit none
	include 'C:\lib\MPI\Include\mpif.h'
	integer :: code, rank
	integer, parameter :: nb_values=1000, tag = 110, dp = kind(1.d0)
  	integer, dimension(MPI_STATUS_SIZE) :: status
  	real(kind=dp), dimension(nb_values)  :: message

	call MPI_INIT(code)
	call MPI_COMM_RANK(MPI_COMM_WORLD,rank,code)


	if (rank == 0) then
     	call random_number(message)
		call MPI_SEND(message,1,MPI_DOUBLE_PRECISION,1,tag,MPI_COMM_WORLD,code)
	elseif (rank==1) then
		call MPI_RECV(message,1,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,status,code)
	print('("Me, process 1, received ",i4," values (last = ", &
            & f4.2,") from process 0.")'), nb_values,message(nb_values) 	
	end if 

	call MPI_FINALIZE(code)

	

end program TP1