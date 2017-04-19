#include <stdio.h>
#include <string.h>
#include <omp.h>
#include <stdlib.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(int argc, char *argv[]) {
	int comm_sz = atoi(argv[2]);
	char* pEnd;
	int N = strtol(argv[1],&pEnd,10);
	int localN = N/comm_sz;//assuming divisible
	double *shared =(double *)malloc((2*comm_sz-2) *sizeof(double));
	
	int end = 120;//end=20N is roughly 1 period
	//int writeoutput = 1;//0 for false
	
	if (N%comm_sz!=0) {
		printf("ERROR N must be divisible by comm_sz. Terminating.");
		exit(1);
	}
	if (localN<2) {
		printf("ERROR N (%d) must be at least twice as large as number of processes (%d). Terminating.",N,comm_sz);
		exit(1);
	}
	
	//double forOutput[N];
	
	#pragma omp parallel num_threads(comm_sz)
	{
	double * f0 = (double *)malloc(localN*sizeof(double));
	double * f1 = (double *)malloc(localN*sizeof(double));
	double * f2 = (double *)malloc(localN*sizeof(double));
	double *temp;
	int my_rank = omp_get_thread_num();
	//printf("my thread %d\n",my_rank);
	double localx0 = 1.0/(N-1)*my_rank*localN;
	//int local_First_Index = my_rank * (N/comm_sz);//use inclusively
	int i;
		#pragma omp for
		for(i=0; i<N; i++){
			double x = localx0 + (double)i*1.0/(N-1);
			f0[i] = sin(M_PI*x);
			f1[i] = sin(M_PI*x);
			//printf("f0 %d %lf\nf1 %d %lf\n",i,f0[i],i,f1[i]);
		}//for
		
		if(my_rank == 0){
			f0[0] = 0;
			f1[0] = 0;
			f2[0] = 0;
		}//if 
		
		if(my_rank == omp_get_num_threads()-1){
			f0[localN-1] = 0;
			f1[localN-1] = 0;
			f2[localN-1] = 0;
		}//if

		double leftNeighbor=0;
		double rightNeighbor=0;
		//int partnerRight = my_rank+1;
		//int partnerLeft = my_rank-1;

		int step = 2;//current step
		while (step<=end) {
			int i;
			for(i=1; i<localN-2; i++){
				f2[i] = 0.01*(f1[i-1]-2*f1[i]+f1[i+1]) + 2*f1[i] - f0[i];
			}//easy part
			
			#pragma omp critical
			{
				if(my_rank != 0){
					*(shared + (my_rank*2-1)) = *f1;
					printf("Left val %lf by thread %d\n",*f1,my_rank);
				}
				if(my_rank != comm_sz-1){
					*(shared + (my_rank*2)) = *(f1 + (localN -1));
					printf("Right val %lf by thread %d\n",*(f1 + (localN -1)),my_rank);
				} 
			}//critical

				//compute left edge of my domain
				if (my_rank!=0) {
					int i=0;
					leftNeighbor = *(shared + (my_rank * 2 - 2));
					f2[i] = 0.01*(leftNeighbor-2*f1[i]+f1[i+1]) + 2*f1[i] - f0[i];
				}
				//compute right edge of my domain
				if (my_rank!=comm_sz-1) {
					int i=localN-1;
					rightNeighbor = *(shared + (my_rank * 2 + 1));
					f2[i] = 0.01*(f1[i-1]-2*f1[i]+rightNeighbor) + 2*f1[i] - f0[i];
				}
		
				if (my_rank==0) {
					//writerow(N,forOutput);
				}
			
			//rearrange pointers for next step
			temp=f0;
			f0=f1;
			f1=f2;
			f2=temp;

			step++;
		}//while
	}//pragma omp parallel num_threads(comm_sz)
	
	return 0;
}
