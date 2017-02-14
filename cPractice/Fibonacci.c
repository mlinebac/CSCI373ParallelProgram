#include <stdio.h>

int fib(int n);

int main(int argc, int *argv []){

printf("fib(%d)=%d\n", 10, fib(10));
return 0;


}
int fib(int n){
 int i;
  if ((n==0) | (n==1)) return n;
  int a = 0;
  int b = 1;
  int temp;
   for (i=0; i<n; i++){
	temp = b;
	b = a+b;
	a = temp;
 }
return a;
}
