! TD MPI
program TD1
	implicit none
	include 'C:\lib\MPI\Include\mpif.h'

	! call MPI_INIT(code)
	! integer :: nb_procs, rank, code
	
	! All the communicators processes
	! MPI_COMM_WORLD

	! Stop excution before the normal ending, not the Fortran instruction stop
	! MPI_ABORT()(comm,error,code)
	! integer, intent(in) :: comm,error
	! integer, intent(out): code

	! Number of processes managed 
	! call MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,code)

	! Obtain the process rank
	! call MPI_COMM_RANK(MPI_COMM_WORLD,rank,code)
	! integer, intent(out) :: rank, code
	! integer, intent(in) :: comm

	! print *, 'I am the process',rank, 'among',nb_procs


	!Block Send MPI_SEND(buf,count,datatype,dest,tag,comm,code)
	!!<type> :: buf
	! integer :: count, datatype 
	! integer :: dest, tag, comm, code
	! integer :: dimension(MPI_STATUS_SIZE) :: status_msg

	!Block Receive MPI_RECEV(buf,count,datatype,source,tag,comm,status_msg,code)
	! <type> :: buff
	! integer :: count,datatype,
	! integer :: source, tag, comm, code
	! integer, dimension(MPI_STATUS_SIZE) :: status_msg

	! integer, dimension(MPI_STATUS_SIZE) :: status_msg
	! integer, parameter :: tag = 100
	! integer :: rank, value, code

	! call MPI_INIT(code)

	! call MPI_COMM_RANK(MPI_COMM_WORLD, rank, code)

	! if (rank == 2) then 
	! 	value = 1000
	! 	call MPI_SEND(value,1,MPI_INTEGER,3,tag,MPI_COMM_WORLD,code)
	! elseif (rank == 3) then
	! 	call MPI_RECV(value,1,MPI_INTEGER,2,tag,MPI_COMM_WORLD,status_msg,code)
	! 	print*, 'I, process 3, I received', value, ' from process 2'
	! end if

	!!! Notes: the rank of the sender and the rank of the tag can be replaced
	!!!!!!!!!! respectively by the MPI_ANY_SOURCE and MPI_ANY_TAG wildcards
	! MPI_PROC_NULL : a dummy process of rank
	! MPI_STATUS_IGNORE is a predefined constant which can be used instead of the status variale


	! MPI_SENDRECV(sendbuf,sendcount,sendtype,dest,sendtag,recvbuf,recvcount
	!,recvtype,source,recvtag,comm,status_msg,code)
	! <type> :: sendbuf, recvbuf
	! integer:: sendcount, recvcount
	! integer :: sendtype, recvtype 
	! integer :: source, dest, sendtage, recvtag, comm, code
	! integer, dimension(MPI_STATUS_SIZE) :: status_msg

	! integer :: rank, value, num_proc, code
	! integer, parameter :: tag = 110

	! call MPI_INIT(code)
	! call MPI_COMM_RANK(MPI_COMM_WORLD, rank, code)

	! num_proc = mod(rank+1,2)
	! ! value = rank+1000
	! call MPI_SENDRECV(rank+1000, 1, MPI_INTEGER,num_proc,tag,&
	! 				value,1,MPI_INTEGER,num_proc,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,code)
	
	! ! call MPI_SEND(rank+1000, 1, MPI_INTEGER, num_proc,tag,MPI_COMM_WORLD,code)
	! ! call MPI_RECV(value, 1, MPI_INTEGER, num_proc, tag, MPI_COMM_WORLD,MPI_STATUS_IGNORE,code)
	! print *, 'I, process', rank, ' I received', value, ' from process', num_proc

	! MPI_SENDRECV_REPLACE(buff, count, datatype,
	!                      dest, sendtag, source, recvtag, comm, status_msg, code)


	! Global synchronization
	! call MPI_BARRIER()

	! Global distribution
	! call MPI_BCAST(buf, count, datatype, root, comm, code)

	! integer :: rank, value, code

	! call MPI_INIT(code)
	! call MPI_COMM_RANK(MPI_COMM_WORLD, rank, code)

	! if (rank==2) value = rank + 1000
	! call MPI_BCAST(value,1,MPI_INTEGER,2, MPI_COMM_WORLD,code)

	! print*, value

	! call MPI_SCATTER(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm, code)

	! integer :: rank, value, code, nb_procs, block_length, i 
	! integer, parameter :: nb_values = 8 
	! real, allocatable, dimension(:) :: values, recvdata

	! call MPI_INIT(code)
	! call MPI_COMM_RANK(MPI_COMM_WORLD, rank, code)
	! call MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,code)

	! block_length =  nb_values/nb_procs
	! allocate(recvdata(block_length))

	! if (rank == 2) then
	! 	allocate(values(nb_values))
	! 	values(:) = (/(1000.+i,i=1,nb_values)/)
	! 	print*, 'I process', rank, ' send my values array', values(1:nb_values)
	! end if

	! call MPI_SCATTER(values,block_length,MPI_REAL,recvdata,block_length,&
	! 				MPI_REAL,2,MPI_COMM_WORLD,code)
	! print*, 'I process', rank, 'received ', recvdata(1:block_length), &
	! 		' of process 2'

	! call MPI_GATHER(Sendbuf,sendcount,sendtype,recvbuf,recvcount,recvtype,root,comm,code)

	! integer, parameter :: nb_values = 8
	! integer :: nb_procs, rank, block_length, i, code
	! real, dimension(nb_values) :: recvdata
	! real, allocatable, dimension(:) :: values

	! call MPI_INIT(code)
	! call MPI_COMM_SIZE(MPI_COMM_WORLD, nb_procs, code)
	! call MPI_COMM_RANK(MPI_COMM_WORLD, rank, code)

	! block_length = nb_values/nb_procs
	! allocate(values(block_length))

	! values(:) = (/(1000.+rank*block_length+i, i = 1, block_length)/)
	! print*, 'I process ', rank, ' sent my values array: ', &
	! 		values(1:block_length)
	! call MPI_GATHER(values,block_length,MPI_REAL,recvdata,block_length,MPI_REAL,2,MPI_COMM_WORLD,code)

	! if (rank == 2) print*,'I process 2 received', recvdata(1:nb_values)



	! call MPI_ALLGATHER(sendbuf,sendcount,sendtype,recvbuf,recvcount,recvtype,comm,code)

	! integer, parameter :: nb_values = 8
	! integer :: nb_procs, rank, block_length, i, code
	! real, dimension(nb_values) :: recvdata
	! real, allocatable, dimension(:) :: values

	! call MPI_INIT(code)
	! call MPI_COMM_SIZE(MPI_COMM_WORLD, nb_procs, code)
	! call MPI_COMM_RANK(MPI_COMM_WORLD, rank, code)

	! block_length = nb_values/nb_procs
	! allocate(values(block_length))

	! values(:) = (/(1000.+rank*block_length+i,i=1,block_length)/)

	! call MPI_ALLGATHER(values,block_length,MPI_REAL,recvdata,block_length,MPI_REAL,MPI_COMM_WORLD,code)
	! print*, 'I process ', rank, 'received', recvdata(1:nb_values)

	! call MPI_GATHERV(sendbuf,sendcount,recvbuf,recvcounts,displs,recvtype,root,comm,code)
	! The i-th process of the comm sends to process root: a message starting at position
	! sendbuf, of sendcount element of type sendtype, and receives at position recvbuff, of
	! recvcounts(i) element of type recvtype with a displacement of displs(i)


	! integer, parameter :: nb_values =10
	! integer :: remainder
	! integer :: nb_procs, rank, block_length, i, code
	! real, dimension(nb_values) :: recvdata
	! real, allocatable, dimension(:) :: values
	! integer, allocatable, dimension(:) :: nb_elements_received, displacement

	! call MPI_INIT(code)
	! call MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,code)
	! call MPI_COMM_RANK(MPI_COMM_WORLD,rank,code)

	! block_length = nb_values/nb_procs
	! remainder = mod(nb_values,nb_procs)
	! if (remainder > rank)	block_length = block_length + 1
	

	! allocate (values(block_length))

	! values(:) = (/(1000.+(rank*(nb_values/nb_procs))+min(rank,remainder)+i,i=1,block_length)/)

	! print*, 'I process ', rank, 'sent my values array ', values(1:block_length)

	! IF (rank == 2 ) THEN
	! 	allocate(nb_elements_received(nb_procs),displacement(nb_procs))
	! 	nb_elements_received(1) = nb_values/nb_procs
	! 	print*, 'nb_elements_received(1)', nb_elements_received(1)
	! 	if (remainder > 0) nb_elements_received(1) = nb_elements_received(1) + 1
	! 	displacement(1) = 0
	! 		DO i = 2, nb_procs
	! 			displacement(i) = displacement(i-1) + nb_elements_received(i-1)
	! 			nb_elements_received(i) = nb_values/nb_procs
	! 			if (remainder > i-1) nb_elements_received(i) = nb_elements_received(i)  + 1
	! 		END DO
	! END IF

	! call MPI_GATHERV(values,block_length,MPI_REAL,recvdata,nb_elements_received,displacement,MPI_REAL,2,MPI_COMM_WORLD,code)
	! if (rank == 2) print*, 'I process 2 received' , recvdata(1:nb_values) 
	

	! call MPI_ALLTOALL(Sendbuf,sendcount,datatype,recvbuf,recvcount,recvtype,comm,code)
	


	! integer, parameter :: nb_values = 8
	! integer :: nb_procs, rank, block_length, i, code
	! real, dimension(nb_values) :: recvdata,values

	! call MPI_INIT(code)
	! call MPI_COMM_RANK(MPI_COMM_WORLD,rank,code)
	! call MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,code)

	! values(:) = (/(1000.+rank*nb_values+i,i=1,nb_values)/)
	! ! block_length =  nb_values/nb_procs
	! block_length = 1
	! print*, 'I process ', rank, ' sent my values array', values
	! call MPI_ALLTOALL(values,block_length,MPI_REAL,recvdata,block_length,MPI_REAL,MPI_COMM_WORLD,code)
	! print*, 'I process ', rank, ' received ', recvdata(1:nb_values)


	! call MPI_SUM() 
	! call MPI_PROD()
	! call MPI_MAX()
	! call MPI_MIN()
	! call MPI_MAXLOC()
	! call MPI_MINLOC()
	! call MPI_LAND()
	! call MPI_LOR()
	! call MPI_LXOR()

	! call MPI_REDUCE()
	! integer :: nb_procs, rank, code, value,sum
	! call MPI_INIT(code)	
	! call MPI_COMM_SIZE(MPI_COMM_WORLD, nb_procs, code)
	! call MPI_COMM_RANK(MPI_COMM_WORLD, rank, code)

	! if (rank==0) then
	! 	value = 1000
	! else
	! 	value = rank
	! end if

	! call MPI_REDUCE(value,sum,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,code)
	! if (rank==0) then
	! 	print*, 'I process 0 have the global value', sum
	! end if


	! call_ALLREDUCE()

	! integer :: nb_procs, rank, code, value, product
	! call MPI_INIT(code)
	! call MPI_COMM_SIZE(MPI_COMM_WORLD, nb_procs, code)
	! call MPI_COMM_RANK(MPI_COMM_WORLD, rank, code)

	! if (rank == 0) then
	! 	value = 10
	! else
	! 	value = rank
	! end if

	! call MPI_ALLREDUCE(value,product,1,MPI_INTEGER,MPI_PROD,MPI_COMM_WORLD,code)

	! print*, 'I process', rank , 'received the value of the global product ', product


	! call MPI_SCAN() subroutine allows making partial reductions by considering, for each
	! , the previsous processes of the communicator and itself
	! call MPI_EXSCAN() is the exclusive version of MPI_SCAN() which is inclusive

	! Personal reduction operations
	! MPI_OP_CREATE()
	! MPI_OP_FREE() 

	! MPI_INPACE keep the result in the same place as the sending buffer

	
	! call MPI_SCAN(sendbuf,recvbuf,count,datatype,op,comm)


	! Communication Modes 
	!! Blocking: MPI_SEND(), MPI_SSEND(), MPI_BSEND(), MPI_RECV()
	!! Non-blocking	: MPI_ISEND(), MPI_ISSEND(), MPI_IBSEND(), MPI_IRECV()

	! A call is blocking if the memory space used for the communication can be reused immediatyly
	! after the exit of the call
	! The data sent can be modified after the call
	! The data receive can be read after the call

	!! Synchronous Sends
	! There can be no communication before the 2 processes are ready to communicate
	! call MPI_SSEND(values, count, msgtype, dest, tag, comm, code)
	! integer :: rank,value,num_proc,code
	! integer,parameter :: tag=110

 ! 	call MPI_INIT(code)
 ! 	call MPI_COMM_RANK( MPI_COMM_WORLD,rank,code)

 ! 	! We run on 2 processes
 ! 	num_proc=mod(rank+1,2)

 ! 	call MPI_SEND(rank+1000,1,MPI_INTEGER,num_proc,tag, MPI_COMM_WORLD,code)
 ! 	call MPI_RECV(value,1, MPI_INTEGER, num_proc,tag, MPI_COMM_WORLD, &
 ! 	MPI_STATUS_IGNORE,code)

	! print*, 'I process ', rank, ' received ', value, ' from process ', num_proc

	
	
	call MPI_FINALIZE(code)
END PROGRAM TD1