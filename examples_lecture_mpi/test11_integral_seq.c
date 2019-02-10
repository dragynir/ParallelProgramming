#include<stdio.h>
#include<math.h>

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

  double a = 0;
  double b = 3.141592653589793238462643383;
  int n = 1000000;

  if (argc>1) sscanf(argv[1],"%d",&n);

  double result = integral(n,a,b);

  printf("N: %10d   Result: %.18lf\n",n,result);
  return 0;
}
