#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(void) {
	
	int test_array[] = {3,4,7,6,13,40,70,60,99,100};//402
	int n = 10;//number of elements
	int m = 4;
	
	int comm_sz;
	int my_rank;
	
	MPI_Init(NULL,NULL);
	MPI_Comm_size(MPI_COMM_WORLD,&comm_sz);
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
	
	//send array data from proc 1 to proc 0
	if (my_rank == 1) {
		//suppose for some reason I might send less than
		//full array
		//MPI_Send(test_array,n,MPI_INT,0,0,MPI_COMM_WORLD);
		MPI_Send(test_array,m,MPI_INT,0,0,MPI_COMM_WORLD);
	} else {
		//but proc 0 doesn't know how many
		//(does know max possible)
		int receive_array[n];
//		MPI_Recv(&receive_array,n,MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//		printf("Contents of receive_array:\n");
//		for (int i=0; i < n; i++) {
//			printf("%d ",receive_array[i]);
//		}
//		printf("\n");
//		MPI_Recv(&receive_array,m,MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//		printf("Contents of receive_array:\n");
//		for (int i=0; i < m; i++) {
//			printf("%d ",receive_array[i]);
//		}
//		printf("\n");

		//determine size of message
		MPI_Status status;
		MPI_Recv(&receive_array,n,MPI_INT, 1, 0, MPI_COMM_WORLD, &status);
		int msg_size;
		MPI_Get_count(&status,MPI_INT,&msg_size);
		printf("message size was: %d\n",msg_size);
		printf("Contents of receive_array:\n");
		for (int i=0; i < msg_size; i++) {
			printf("%d ",receive_array[i]);
		}
		printf("\n");

	}

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
