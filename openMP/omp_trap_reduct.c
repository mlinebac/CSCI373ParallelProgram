#include <stdio.h>
#include <string.h>
#include <omp.h>

const int MAX_STRING = 100;

double f(double x) {
  return x*x;
}

int main(int argz, char **argv) {
  
 
 
  double a=1.0;
  double b=3.0;
  int n = 1000000000;
  double h = (b-a)/(double)n;
  double sum = (f(a)+f(b))/2.0;
  
  
 #pragma omp parallel for reduction(+: sum)
  for (int i=0; i<=n-1; ++i) {
    double x_i = a + i*h;
    sum += f(x_i);
  }
  
  double solution = h*sum;
  
  printf("solution %f\n", solution);
  return 0;
}
