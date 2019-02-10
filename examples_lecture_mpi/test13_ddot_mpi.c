#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<mpi.h>

double init_x(int i,int n) { return (double)i; }
double init_y(int i,int n) { return (double)(n-i); }

double ddot(int n,double *x,double *y) {
  double s = 0.0;
  for (int i=0;i<n;i++) s += x[i] * y[i];
  return s;
}

int main(int argc,char *argv[]) {
  int size,rank;
  int n = 1000;
  if (argc>1) sscanf(argv[1],"%d",&n);

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);

  int k1 = n*rank/size;
  int k2 = n*(rank+1)/size;
  int local_n = k2 - k1;

  double *x = (double *)malloc(local_n*sizeof(double));
  double *y = (double *)malloc(local_n*sizeof(double));

  for (int i=0;i<local_n;i++) x[i] = init_x(k1+i,n);
  for (int i=0;i<local_n;i++) y[i] = init_y(k1+i,n);

  double r = ddot(local_n,x,y);

  double result;
  MPI_Reduce(&r,&result,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

  if (rank == 0) printf("N: %10d   NP: %2d   Result: %.16lf\n",n,size,result);
  free(x); free(y);
  MPI_Finalize();
  return 0;
}
