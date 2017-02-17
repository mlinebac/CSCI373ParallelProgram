#include <stdio.h>
#include <stdlib.h>

int main(void){

int n = 3;
char x[] = "XML";

int i;

printf("%s\n", x);

for(i=n-1; i>=0; i--){
  printf("%c",x[i]);
}
printf("\n");

return 0;
}
