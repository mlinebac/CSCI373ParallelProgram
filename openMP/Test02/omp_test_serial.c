#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

double initialCondition(double x, double y);
void writeheader(int N, int end);
void writerow(int N, double rawdata[]);
void printArray(double **array, int local_N, int N);


int main(int argc, char *argv[]) {
	int comm_sz = 1;
	int my_rank = 0;
	int N = 10;
	int local_N;
	int end = 2;//end=20N is roughly 1 period
	//int writeoutput = 1;//0 for false
	double x,y,local_x, local_y;
	int i,j;
	
	//comm_sz = omp_get_num_threads();
	//my_rank = omp_get_thread_num();
	local_N = N/comm_sz;//assuming divisible
	
	if (N%comm_sz!=0) {
		printf("ERROR N must be divisible by thread_count. Terminating.");
		exit(1);
	}
	if (local_N<2) {
		printf("ERROR N (%d) must be at least twice as large as number of processes (%d). Terminating.",N,comm_sz);
		exit(1);
	}
	
	double** f0 = (double**)malloc(sizeof(double*)*N);
	double** f1 = (double**)malloc(sizeof(double*)*N);
	double** fend = (double**)malloc(sizeof(double*)*N);
	
	
	f0[0] = (double*)malloc(sizeof(double) * local_N * N);
	f1[0] = (double*)malloc(sizeof(double) * local_N * N);
	fend[0] = (double*)malloc(sizeof(double) * local_N * N);
	
	for (i=1; i<local_N; i++){
		f0[i] = (*f0 + local_N * i);
		f1[i] = (*f1 + local_N * i);
		fend[i] = (*fend + local_N * i);
	}
	
	double **originalf0 = f0;
	double **originalf1 = f1;
	double **originalfend = fend;
	
		
	local_x = 1.0/(N-1)*my_rank*local_N;
	local_y = 1.0/(N-1)*my_rank*local_N;
	
	for	(i=0; i < local_N; i++){
		for (j=0; j < N; j++){
			x = local_x + (double)i*1.0/(N-1);
			y = local_y + (double)j*1.0/(N-1);
			f0[i][j] = initialCondition(x,y);
			f1[i][j] = initialCondition(x,y);
		}
	}
	//printArray(f1,local_N, N);
	
	f0 = originalf0;
	f1 = originalf1;
	fend = originalfend;
	
	
	free(f0);
	free(f1);
	free(fend);
	
	
	return 0;
}








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
	fp = fopen("outfile.pgm", "w");
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

void writerow(int N, double rawdata[]) {
	FILE *fp;
	fp = fopen("outfile.pgm", "a");
	if (fp == NULL) {
		printf("sorry can't open outfile.pgm. Terminating.\n");
		exit(1);
	}
	else {
		for (int i=0; i<N; i++) {
			int val = rawdata[i]*127+127;
			fprintf(fp,"%d ",val);
		}
		fprintf(fp,"\n");
		fclose(fp);
	}
}
void printArray(double **array, int local_N, int N){
	int i,j;
	for (i=0; i<local_N; i++){
		for(j=0; j<N; j++){
			printf("%f  ", array[i][j]);
		}
			printf("\n");
	}
	printf("\n");
}
