#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(void) {
	
	
	int comm_sz;
	int my_rank;
	
	MPI_Init(NULL,NULL);
	MPI_Comm_size(MPI_COMM_WORLD,&comm_sz);
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
	
	//try MPI_Bcast
	
	int secret_value = 0;
	
	if (my_rank == 0) {
		secret_value = 33;
	}
	
	printf("before MPI_Bcast, %d has %d\n",my_rank,secret_value);
	
	MPI_Bcast(&secret_value,1,MPI_INT,0,MPI_COMM_WORLD);
	
	printf("after MPI_Bcast, %d has %d\n",my_rank,secret_value);

	MPI_Finalize();
	
	return 0;
}

//double array_sum(int * array, int start, int end) {
long long array_sum(int * array, int start, int end) {
	//double sum=0;
	long long sum=0;
	for (int i=start; i<end; ++i) {
		sum+=array[i];
	}
	return sum;
}
