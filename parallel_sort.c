#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>


int compare_ints(const void * a, const void * b){
	return ( *(int*)a - *(int*)b );
}

void merge_low(int my_keys[], int recv_keys[], int temp_keys[], int n);

void merge_high(int my_keys[], int recv_keys[], int temp_keys[], int n);

void odd_even_sort(int *array, int array_size, int my_rank, int comm_sz);

int compute_partner(int phase, int my_rank, int comm_sz);


int main(void){
	int my_rank,comm_sz;
	MPI_Init(NULL,NULL);
	MPI_Comm_size(MPI_COMM_WORLD,&comm_sz);
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
	
	int arr_len = 10;
	int *arr1 = (int*) malloc(arr_len*sizeof(int));
	for(int i = 0; i < arr_len; i++){
		arr1[i] = rand()%arr_len;
	}
	if(my_rank==0){
		for(int i = 0; i < arr_len; i++){
			printf("%d, ",arr1[i]);
		}
		printf("\n");
	}

	odd_even_sort(arr1,arr_len,my_rank,comm_sz);

	if(my_rank==0){
		for(int i = 0; i < arr_len; i++){
			printf("%d, ",arr1[i]);
		}
		printf("\n");
	}

	MPI_Finalize();
	return 0;
}


void odd_even_sort(int *array, int array_size, int my_rank, int comm_sz){
	int std_array_size = array_size/comm_sz;
	int my_array_size = std_array_size;
	if(my_rank == comm_sz-1 && array_size%comm_sz!=0){
		my_array_size += array_size%comm_sz;
	}
	int *my_array = &array[(my_rank*std_array_size)];
	int partner, partner_array[my_array_size], temp_array[my_array_size];
	qsort(my_array,my_array_size,sizeof(int),compare_ints);
	for (int phase = 0; phase <= comm_sz; phase++){
		partner = compute_partner(phase,my_rank,comm_sz);
		/*printf("Rank: %d Phase: %d\n",my_rank,phase);
		for (int i = 0; i < my_array_size; i++){
			printf("%d, ",my_array[i]);
		}
		printf("\n");*/
		if (partner != MPI_PROC_NULL){
			if (my_rank < partner){
				MPI_Send(my_array,std_array_size,MPI_INT,partner,0,MPI_COMM_WORLD);
				MPI_Recv(partner_array,std_array_size,MPI_INT,partner,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				merge_low(my_array,partner_array,temp_array,my_array_size);
			}else{
				MPI_Recv(partner_array,std_array_size,MPI_INT,partner,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				MPI_Send(my_array,std_array_size,MPI_INT,partner,0,MPI_COMM_WORLD);
				merge_high(my_array,partner_array,temp_array,my_array_size);

			}
		}
	}
}


int compute_partner(int phase, int my_rank, int comm_sz){
	int partner = MPI_PROC_NULL;
	if(my_rank%2 != phase%2){
		if(my_rank>0){
			partner = my_rank - 1;
		}
	}else{
		if(my_rank < comm_sz - 1){
			partner = my_rank + 1;
		}
	}
	return partner;
}

void merge_low(int my_keys[], int recv_keys[], int temp_keys[], int n){
	int m_i, r_i, t_i;
	m_i = r_i = 0;
	for(t_i = 0; t_i < n; t_i++){
		temp_keys[t_i] = (my_keys[m_i]<=recv_keys[r_i])? my_keys[m_i++] : recv_keys[r_i++];
	}
	for(m_i = 0; m_i < n; m_i++){
		my_keys[m_i] = temp_keys[m_i];
	}
}

void merge_high(int my_keys[], int recv_keys[], int temp_keys[], int n){
	int m_i, r_i, t_i;
	m_i = r_i = n-1;
	for(t_i = n-1; t_i >= 0; t_i--){
		temp_keys[t_i] = (my_keys[m_i]>=recv_keys[r_i])? my_keys[m_i--] : recv_keys[r_i--];
	}
	for(m_i = 0; m_i < n; m_i++){
		my_keys[m_i] = temp_keys[m_i];
	}
}
