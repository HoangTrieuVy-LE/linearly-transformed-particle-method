program TP1 
	implicit none
	include 'C:\lib\MPI\Include\mpif.h'

||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

	! integer :: code, rank

	! call MPI_INIT(code)
	! call MPI_COMM_RANK(MPI_COMM_WORLD,rank,code)


	! if (mod(rank,2)/=0) then
	! 	print*, 'I am the odd-ranked process, my rank is ', rank
	! else 
	! 	print*, 'I am the even-ranked process, my rank is ', rank
	! end if

||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


	call MPI_FINALIZE(code)

end program TP1