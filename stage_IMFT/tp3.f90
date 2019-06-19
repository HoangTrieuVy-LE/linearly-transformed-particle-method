program transposee_dune_matrice
implicit none

include 'mpif.h'

integer :: rank,code,type_line,type_transpose,size_real,i,j
integer, parameter : nb_ligne=5,nb_collone=4,tag =100
real, dimension(nb_ligne,nb_collone): A
real, dimension(nb_collone,nb_ligne): AT
integer(kind=MPI_ADRESS_KIND) :: sizedisplacement
integer,dimension(MPI_STATUS_SIZE) :: status

call MPI_INIT(code)
call MPI_COMM_RANK(MPI_COMM_WOWRLD,rank, code)
call MPI_TYPE_SIZE(MPI_REAL,size_real,code)
call MPI_TYPE_VECTOR(nb_collone,1,sizedisplacement,type_line,type_transpose,code)
call


call MPI_INIT(code)

end program transpose_dune_matrice
