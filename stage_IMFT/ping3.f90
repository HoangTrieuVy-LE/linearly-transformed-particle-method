program ping3

! include 'C:\lib\MPI\Include\mpif.h'

include 'mpif.h'

integer :: code, rank
integer, dimension(MPI_STATUS_SIZE) :: status
integer, parameter :: nb_values_max=70000000, nb_boucle = 1000, dp = kind(1.d0),nb_tests=9
integer, dimension(nb_tests)                :: nb_values
real(kind = dp), dimension(0:nb_values_max-1) :: message
real(kind = dp) :: time_begin, time_end

nb_values = (/0,1,10,100,1000,10000,100000,1000000,5000000/)
call MPI_INIT(code)

call MPI_COMM_RANK(MPI_COMM_WORLD, rank, code)

do i= 1,9
	if (rank==0) then
		time_begin = MPI_WTIME()
		call random_number(message)
		call MPI_SEND(message,nb_values(i),MPI_DOUBLE_PRECISION,1,tag,MPI_COMM_WORLD,code)
		call MPI_RECV(message,nb_values(i),MPI_DOUBLE_PRECISION,1,tag,MPI_COMM_WORLD,status,code)
		time_end = MPI_WTIME()
		call print_function()
	else if (rank == 1) then
		call MPI_RECV(message,nb_values(i),MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,status,code)
		call MPI_SEND(message,nb_values(i),MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,code)
	end if
end do

call MPI_FINALIZE(code)

contains
	subroutine print_function
	if(nb_values(i)/=0) then
		print('("Me, process 0, sent and received ",&
i8,"values (last = ",f8.6,")from process 1 in ",f8.6,"second,bandwidth",f7.2,"Mo/s")'),&
nb_values(i),message(nb_values(i)-1),time_end-time_begin,real(2*nb_values(i)*8)/1000000./(time_end-time_begin)


	else
       		print ('("Me, process 0, sent and received ",i8, &
            & " values in ",f8.6," seconds, bandwidth ",f7.2, &
            & " Mo/s.")'), &
            nb_values(i),time_end-time_begin, &
            real(2*nb_values(i)*8)/1000000./(time_end-time_begin)
    	end if
	end subroutine print_function
end program ping3
