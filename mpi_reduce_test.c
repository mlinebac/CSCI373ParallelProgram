#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(void) {
	
	int test_array[] = {3,4,7,6,13,40,70,60,99,100};//402
	int n = 10;//number of elements
	int comm_sz;
	int my_rank;
	
	MPI_Init(NULL,NULL);
	MPI_Comm_size(MPI_COMM_WORLD,&comm_sz);
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
	
	//try MPI_Reduce
	if (my_rank != 0) {
		MPI_Reduce(&test_array[my_rank],NULL,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
	}
	else {
		int sum;
		MPI_Reduce(&test_array[my_rank],&sum,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
		printf("sum is %d\n",sum);
	}

	MPI_Finalize();
	
	return 0;
}

