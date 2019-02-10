#include<stdio.h>
#include<mpi.h>

#define N 10
int x[N];

int main(int argc,char *argv[]) {
  int size,rank;
  MPI_Status st;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  for (int i=0;i<N;i++) x[i] = rank;

  if (size >= 3) {
    if (rank == 2) MPI_Send(x,N,MPI_INT,0,12345,MPI_COMM_WORLD);
    if (rank == 0) MPI_Recv(x,N,MPI_INT,2,12345,MPI_COMM_WORLD,&st);
  }

  printf("Hello! I am %d of %d, Data: %d\n",rank,size,x[0]);
  MPI_Finalize();
  return 0;
}
