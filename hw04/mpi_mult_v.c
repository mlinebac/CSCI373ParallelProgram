//Matt Lineback
//CSCI 373 Parallel Programming HW04
//MPI program that implements multiplication of a vector
//by a scalar and dot product.

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>


int main(void){

int n,local_n;
int i;
int scalar;
int *V1, *V2;
int *local_V1, *local_V2;
int local_sum;
int sum;
int comm_sz;
int my_rank;

MPI_Init(NULL,NULL);
MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

 if(my_rank == 0){
	printf("enter vector size \n");
	scanf("%d", &n);
	local_n = n/comm_sz;
	}

MPI_Bcast (&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
MPI_Bcast (&local_n, 1, MPI_INT, 0, MPI_COMM_WORLD);
		
 if(my_rank == 0){
	printf("enter scalar \n");
	scanf("%d", &scalar);
MPI_Bcast (&scalar, 1, MPI_INT, 0, MPI_COMM_WORLD);
	}

MPI_Bcast (&scalar, 1, MPI_INT, 0, MPI_COMM_WORLD);

 V1 = (int *) malloc (n * sizeof(int));
 V2 = (int *) malloc (n * sizeof(int));
 local_V1 = (int *) malloc (local_n * sizeof(int));
 local_V2 = (int *) malloc (local_n * sizeof(int));

	
 if (my_rank == 0){
   printf("enter the first vector\n");
	for (i=0; i<n; i++){
	 scanf("%d", &V1[i]);
	}
   printf("enter the second vector\n");
	for (i=0; i<n; i++){
	 scanf("%d", &V2[i]);
	}

MPI_Scatter(V1, local_n, MPI_INT, local_V1, local_n, MPI_INT, 0, MPI_COMM_WORLD);

MPI_Scatter(V2, local_n, MPI_INT, local_V2, local_n, MPI_INT, 0, MPI_COMM_WORLD);
	
  free(V1);
  free(V2);

 }else{

MPI_Scatter(V1, local_n, MPI_INT, local_V1, local_n, MPI_INT, 0, MPI_COMM_WORLD);

MPI_Scatter(V2, local_n, MPI_INT, local_V2, local_n, MPI_INT, 0, MPI_COMM_WORLD);
	
 }

 for (i=0; i<local_n; i++){
  local_V1[i] = local_V1[i] * scalar;
  local_V1[i] = local_V1[i] * local_V2[i];
  local_sum = local_sum + local_V1[i];
 }

MPI_Reduce (&local_sum, &sum, 1, MPI_INT, MPI_SUM, 0,
		 MPI_COMM_WORLD);

 if(my_rank == 0){
  printf("the result is: %d\n", sum);

  free(local_V1);
  free(local_V2);

 }

	
MPI_Finalize();
return 0;
}

