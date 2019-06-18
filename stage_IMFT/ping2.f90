program TP2
	implicit none
	include 'C:\lib\MPI\Include\mpif.h'

	integer :: code, rank
	integer, parameter :: nb_values=1000, tag = 110, dp = kind(1.d0)
  	integer, dimension(MPI_STATUS_SIZE) :: status
  	real(kind=dp), dimension(nb_values)  :: message

	call MPI_INIT(code)
	call MPI_COMM_RANK(MPI_COMM_WORLD, rank, code)
	if (rank == 0) then
		call random_number(message)
	     time_begin=MPI_WTIME()
		call MPI_SEND(message,nb_values,MPI_DOUBLE_PRECISION,1&
					,MPI_COMM_WORLD,tag,code)
		call MPI_RECV(message,nb_values,MPI_DOUBLE_PRECISION,1&
					,MPI_COMM_WORLD,tag,status,code)
	    time_end=MPI_WTIME()
		print ('("Me, process 0, sent and received ",i5, &
         & " values (last = ",f4.2,") from process 1", &
         & " in ",f8.6," seconds.")'), &
           nb_values,values(nb_values),time_end-time_begin
	elseif ( rank == 1 ) then
		call MPI_RECV(message,nb_values,MPI_DOUBLE_PRECISION,0&
					,MPI_COMM_WORLD,tag,status,code)
		call MPI_SEND(message,nb_values,MPI_DOUBLE_PRECISION,0,&
					MPI_COMM_WORLD,tag,code)
	end if 



	call MPI_FINALIZE(code)

end program TP2