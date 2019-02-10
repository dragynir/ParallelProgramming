#include<stdio.h>
#include<math.h>
#include<mpi.h>

double f(double x) { return exp(x)*sin(x); }

double integral(int n,double a,double b) {
  double h = (b-a)/n;
  double s = 0;
  double f1 = f(a);
  for (int i=0;i<n;i++) {
    double x = a + (i+1)*h;
    double f2 = f(x);
    s += 0.5*(f1 + f2);
    f1 = f2;
  }
  return h*s;
}

int main(int argc,char *argv[]) {
  int size,rank;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  double a = 0;
  double b = 3.141592653589793238462643383;
  int n = 1000000;

  if (rank == 0)
    if (argc>1) sscanf(argv[1],"%d",&n);

  MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);

  int k1 = n*rank/size;
  int k2 = n*(rank+1)/size;
  int local_n = k2 - k1;
  double h = (b-a)/n;
  double local_a = a + k1*h;
  double local_b = a + k2*h;
  double r = integral(local_n,local_a,local_b);

  double result;
  MPI_Reduce(&r,&result,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

  if (rank == 0) printf("N: %10d   NP: %2d   Result: %.16lf\n",n,size,result);
  MPI_Finalize();
  return 0;
}
