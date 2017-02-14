//Matt Lineback
//CSCI 373 Parallel Programming HW04
//Program Assignment 3.2 monte carlo method and pi estimation

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

int main(void){

long long int num_in_circle;
long long int local_num_in;
int num_of_toss;
int local_num_toss;
int i;
double dist_sqr;
double x,y;
double pi_est;
int comm_sz;
int my_rank;

MPI_Init(NULL, NULL);
MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

 if (my_rank == 0){
	printf("enter number of tosses: \n");
	scanf("%d", &num_of_toss);
	local_num_toss = num_of_toss/comm_sz;
 }

MPI_Bcast (&num_of_toss, 1, MPI_INT, 0, MPI_COMM_WORLD);
MPI_Bcast (&local_num_toss, 1, MPI_INT, 0, MPI_COMM_WORLD);

	for (i=0; i<local_num_toss; i++){
		x = 2*(double)rand()/(RAND_MAX) -1;
		y = 2*(double)rand()/(RAND_MAX) -1; 
		dist_sqr = x*x + y*y;
		if (dist_sqr <= 1){
			local_num_in++;
		}
	}

MPI_Allreduce(&local_num_in, &num_in_circle, 1, MPI_INT, MPI_SUM, 			MPI_COMM_WORLD);
	
pi_est = 4*num_in_circle/((double)num_of_toss);	
		
 if (my_rank == 0){
	printf("number of tosses: %d\n", num_of_toss);
	printf("estimate of Pi: %.25f\n", pi_est);
 }

MPI_Finalize();
return 0;
}


