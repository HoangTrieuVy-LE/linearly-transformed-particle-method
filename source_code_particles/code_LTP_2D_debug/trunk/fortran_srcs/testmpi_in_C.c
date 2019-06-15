#include <C:\lib\MPI\Include\mpi.h>
#include <stdio.h>
#include <unistd.h>

int main(int argc, char **argv)
{
  int rank;
  

  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  

  printf("Hello world!  I am process number: %d on host \n", rank);

  MPI_Finalize();

  return 0;
}