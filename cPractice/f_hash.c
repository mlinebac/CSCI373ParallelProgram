#include <stdio.h>
#include <stdlib.h>

int main(void){

int f[7] = {6, 1, 1, 5, 1, 1, 1};
int n = 7;
char hash[7];
int i,j;

for (i=0; i<n; i++){
	printf("\n");
for (j=0; j<f[i]; j++)
	printf("#");
}
printf("\n");
return 0;
}
