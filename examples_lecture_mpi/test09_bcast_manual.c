#include<stdio.h>
#include<mpi.h>

int main(int argc,char *argv[]) {
  int size,rank;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  int n = 0;

  if (rank==0) {
    printf("Enter the parameter: ");
    scanf("%d",&n);
  }

  if (size>1) {
    if (rank==0) {
      for (int r=1;r<size;r++)
        MPI_Send(&n,1,MPI_INT,r,r,MPI_COMM_WORLD);
    }
    else
      MPI_Recv(&n,1,MPI_INT,0,rank,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  }

  printf("%2d : %d\n",rank,n);

  MPI_Finalize();
  return 0;
}
