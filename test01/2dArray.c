#include<stdio.h>
#include<stdlib.h>
  
int main()
{
    int r=3, c=4;
    int **arr;
    int count = 0,i,j;
  
    arr  = (int **)malloc(sizeof(int *) * r);
    arr[0] = (int *)malloc(sizeof(int) * c * r);
 
    for(i = 0; i < r; i++)
        //arr[i] = (*arr + c * i);
		arr[i] = (int *)malloc(sizeof(int) * c * i);
    for (i = 0; i < r; i++)
        for (j = 0; j < c; j++)
            arr[i][j] = ++count;  // OR *(*(arr+i)+j) = ++count
  
   // for (i = 0; i <  1; i++){
        for (j = 0; j < c; j++){
            printf("%d ", arr[1][j]);
		}
		printf("\n");
	//}
    return 0;
}
