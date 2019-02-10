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

  MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);

  printf("%2d : %d\n",rank,n);

  MPI_Finalize();
  return 0;
}
