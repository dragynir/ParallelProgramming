#include<stdio.h>
#include<mpi.h>

#define N 100000
int x[N],y[N];

int main(int argc,char *argv[]) {
  int size,rank;
  MPI_Request req[2];
  MPI_Status st[2];
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  for (int i=0;i<N;i++) { x[i] = rank; y[i] = -1; }

  if (size > 1) {
    MPI_Status status;

    MPI_Isend(x,N,MPI_INT, (rank-1+size)%size ,12345,MPI_COMM_WORLD,&req[0]);
    MPI_Irecv(y,N,MPI_INT, (rank+1     )%size ,12345,MPI_COMM_WORLD,&req[1]);
    MPI_Waitall(2,req,st);
  }

  printf("%2d: %d\n",rank,y[0]);
  MPI_Finalize();
  return 0;
}
