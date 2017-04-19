/*
 * Author: Matt Lineback
 * CSCI 373 Parallel Programming
 * 4/12/2017
 * Test 02 openMP
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

double initialCondition(double x, double y);
void writeheader(int N, int end);
void writerow(double **rawdata, int N);
void printArray(double **array, int local_N);

int main(int argc, char *argv[]) {
	int N = 8;
	int local_N;
	int comm_sz = 2;//atoi(argv[1]);
	int end = 8;
	int writeoutput = 0;
	double x,y,local_x,local_y;
	int i,j;
	local_N = N/comm_sz;
	double output[N][N];
	double **shared = (double**)malloc(N * sizeof(double*));
	shared[0] = (double*)malloc(N * N * sizeof(double*));
	
	for (i=0; i<N; i++){
		shared[i] = (*shared + N * i);
	}
	
	if (N%comm_sz!=0) {
		printf("ERROR N must be divisible by thread_count. Terminating.");
		exit(1);
	}
	if (local_N<2) {
		printf("ERROR N (%d) must be at least twice as large as number of processes (%d). Terminating.",N,comm_sz);
		exit(1);
	}
	
	double **f0 = (double**)malloc(N * sizeof(double*));
	double **f1 = (double**)malloc(N * sizeof(double*));
	double **fend = (double**)malloc(N * sizeof(double*));
	double **temp;// = (double**)malloc(sizeof(double*)*N);
	
	f0[0] = (double*)malloc(N * local_N * sizeof(double));
	f1[0] = (double*)malloc(N * local_N * sizeof(double));
	fend[0] = (double*)malloc(N * local_N * sizeof(double));

	for (i=0; i<N; i++){
	  f0[i] = (*f0 + local_N * i);
      f1[i] = (*f1 + local_N * i);
      fend[i] = (*fend + local_N * i);
	}
	
	double **originalf0 = f0;
	double **originalf1 = f1;
	double **originalfend = fend;
	
	
	int my_rank = omp_get_thread_num();
	local_x = 1.0/(N-1)*my_rank*local_N;
	local_y = 1.0/(N-1)*my_rank*local_N;
	
	
	for	(i=1; i < N-1; i++){
		for (j=1; j < N-1; j++){
			x = local_x + (double)i*1.0/(N-1);
			y = local_y + (double)j*1.0/(N-1);
			f0[i][j] = initialCondition(x,y);
			f1[i][j] = initialCondition(x,y);
		}
	}
	
	if (writeoutput){
		writeheader(N,end);
	}
	
	double *leftNeighbor = malloc(N * sizeof(double));
	double *rightNeighbor = malloc(N * sizeof(double));
	#pragma omp parallel num_threads(comm_sz)
	{
	int step = 2;
	  while (step<=end) {
		int i,j;
		#pragma omp for
	    for(i=1; i<N-1; i++){
		  for(j=1; j<N-1; j++){
		  fend[i][j] = 0.01*(f1[i-1][j]+f1[i+1][j]+f1[i][j-1]+f1[i][j+1]-4*f1[i][j])+2*f1[i][j]-f0[i][j];
		  }
		}	
		#pragma omp master
		{
		for(i=0; i<N; i++){
		 for(j=0; j<N; j++){
		   output[i][j] = fend[i][j];
			}
		  }
		}
		#pragma omp critical
		{
		if(my_rank != 0){
		  for(i = 0; i < N; i++){
			for(j = local_N; j<N-1; j++){
				shared[i][j] = f1[i][j];
			}
		  }
		}
		if(my_rank != comm_sz-1){
		  for(i = 0; i < N; i++){
			for(j=0; j < local_N; j++){
				shared[i][j] = f1[i][j];
			}
		  } 
		}
	}//end critical
		//compute left edge of my domain
		if (my_rank!=0) {
		  for(i=0; i<N; i++){
	        for(j=local_N; j<N-1; j++){
			 leftNeighbor[i] = shared[i][my_rank * 2 + 1];
			  fend[i][j] = 0.01*(f1[i-1][j]+f1[i+1][j]+leftNeighbor[i]+f1[i][j+1]-4*f1[i][j])+2*f1[i][j]-f0[i][j];
			}
		  }
		}
		//compute right edge of my domain
		if (my_rank!=comm_sz-1) {
		  for(i=0; i<N; i++){
		    for(j=0; j<local_N; j++){
			  rightNeighbor[i] = shared[i][my_rank*2-2];
			  fend[i][j] = 0.01*(f1[i-1][j]+f1[i+1][j]+f1[i][j-1]+rightNeighbor[i]-4*f1[i][j])+2*f1[i][j]-f0[i][j];
			 }
		   }
		}
	//#pragma omp single
	{
		//rearrange pointers for next step
		temp=f0;
		f0=f1;
		f1=fend;
		fend=temp;
	}
		step++;
			
		}//while
	
	f0 = originalf0;
	f1 = originalf1;
	fend = originalfend;

	free(f0);
	free(f1);
	free(fend);
}//end parallel block
	return 0;
}//end main

double initialCondition(double x, double y) {
	//double sigma=0.01;//tight point
	double sigma=0.1;//wider point
	double mu=0.5;//center
	double max = (1.0/(2.0*M_PI*sigma*sigma))*exp(-0.5*( ((0.5-mu)/sigma)*((0.5-mu)/sigma) +  ((0.5-mu)/sigma)*((0.5-mu)/sigma)   ));
	double result = (1.0/(2.0*M_PI*sigma*sigma))*exp(-0.5*( ((x-mu)/sigma)*((x-mu)/sigma) +  ((y-mu)/sigma)*((y-mu)/sigma)   ))/max;
	return result;
}
void writeheader(int N, int end) {
	FILE *fp;
	char outputfile[20];
	int i;
	for (i=0; i<N; i++){
		sprintf(outputfile,"output%04d.pgm", i);
		fp = fopen(outputfile, "w");
		
	if (fp == NULL) {
		printf("sorry can't open outfile.pgm. Terminating.\n");
		exit(1);
	}
	else {
		// print a table header
		fprintf(fp, "%s\n%d %d\n%s\n", "P2", N, end, "255");
		fclose(fp);
	}
  }
}

void writerow(double **rawdata, int N) {
	FILE *fp;
	char outputfile[20];
	int i, j, k;
	for (k=0; k<N; k++){
	sprintf(outputfile,"output%04d.pgm",k);
	fp = fopen(outputfile, "a");

	if (fp == NULL) {
		printf("sorry can't open outfile.pgm. Terminating.\n");
		exit(1);
	}
	else {
		for(i=0; i<N; i++){
		  for (j=0; j<N; j++){
			int val = rawdata[i][j]*254+127;
			fprintf(fp,"%d ",val);
		  }
		  fprintf(fp, "\n");
		  //fclose(fp);
		}
		fclose(fp);
	  }
	  
	}
  }

void printArray(double **array, int local_N){
	int i,j;
	for (i=0; i<local_N; i++){
		for(j=0; j<local_N; j++){
			printf("%f  ", array[i][j]);
		}
			printf("\n");
	}
	printf("\n");
}

