#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double init_x(int i,int n) { return (double)i; }
double init_y(int i,int n) { return (double)(n-i); }

double ddot(int n,double *x,double *y) {
  double s = 0.0;
  for (int i=0;i<n;i++) s += x[i] * y[i];
  return s;
}

int main(int argc,char *argv[]) {

  int n = 1000;
  if (argc>1) sscanf(argv[1],"%d",&n);

  double *x = (double *)malloc(n*sizeof(double));
  double *y = (double *)malloc(n*sizeof(double));

  for (int i=0;i<n;i++) x[i] = init_x(i,n);
  for (int i=0;i<n;i++) y[i] = init_y(i,n);

  double result = ddot(n,x,y);

  printf("N: %10d   Result: %.18lf\n",n,result);
  free(x); free(y);
  return 0;
}
