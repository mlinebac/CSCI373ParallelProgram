#include <stdio.h>

int isprime(int n);

int main(int argc, char *argv[]){
 isprime(12);
 isprime(17);


return 0;

}
 int isprime(int n){
  if(n<=1) return 0;
   int i;
   for(i=2; (i*i)<=n; i++){
    if (n % i == 0) {
	printf("%d is even\n",n);
	return 0;
	}
}
printf("%d is odd\n",n);
return 1; 

}
