#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>
//#ifndef M_PI
//#define M_PI 3.14159265359

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

double initialCondition(double x, double y) {
	//double sigma=0.01;//tight point
	double sigma=0.1;//wider point
	double mu=0.5;//center
	double max = (1.0/(2.0*M_PI*sigma*sigma))*exp(-0.5*( ((0.5-mu)/sigma)*((0.5-mu)/sigma) +  ((0.5-mu)/sigma)*((0.5-mu)/sigma)   ));
	double result = (1.0/(2.0*M_PI*sigma*sigma))*exp(-0.5*( ((x-mu)/sigma)*((x-mu)/sigma) +  ((y-mu)/sigma)*((y-mu)/sigma)   ))/max;
	return result;
}

void fInitOuterBounds(double **array, int localN);

int main(int argc, char *argv[]) {
	int comm_sz;
	int my_rank;
	int N = 12;
	int localN;
	int end = 120;//end=20N is roughly 1 period
	int writeoutput = 1;//0 for false
	
	
	MPI_Init(NULL,NULL);
	MPI_Comm_size(MPI_COMM_WORLD,&comm_sz);
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
	
	localN = N/comm_sz;//assuming divisible
	
	if (N%comm_sz!=0) {
		printf("ERROR N must be divisible by comm_sz. Terminating.");
		exit(1);
	}
	if (localN<2) {
		printf("ERROR N (%d) must be at least twice as large as number of processes (%d). Terminating.",N,comm_sz);
		exit(1);
	}
	
	double **f0 = malloc(N * sizeof *f0);
	double **f1 = malloc(N * sizeof *f1);
	double **fend = malloc(N * sizeof *fend);
	
	int i,j;
	 for(i=0; i<N; i++){
		 f0[i] = malloc(N*sizeof *f0[i]);
		 f1[i] = malloc(N*sizeof *f1[i]);
		 fend[i] = malloc(N*sizeof *fend[i]);
	}
	
	double localx = 1.0/(N-1)*my_rank*localN;
	double localy = 1.0/(N-1)*my_rank*localN;
	
	double **originalf0 = f0;
	double **originalf1 = f1;
	double **originalfend = fend;
	
	double x;
	double y; 
	
	for (int i=1; i<localN-1; i++){
		for (int j=1; j<=localN-2; j++){
			x = localx + (double)i*1.0/(N-1);
			y = localy + (double)j*1.0/(N-1);
			f0[i][j] = initialCondition(x,y);
			f1[i][j] = initialCondition(x,y);
	}
}	
	if (my_rank == 0){
		fInitOuterBounds(f0,N);
		fInitOuterBounds(f1,N);
		fInitOuterBounds(fend,N);
	
	}
	
	if(my_rank == comm_sz-1){
		fInitOuterBounds(f0,localN);
		fInitOuterBounds(f1,localN);
		fInitOuterBounds(fend,localN);
	}
	
	
	for (int i=0; i<localN; i++){
		for(int j=0; j<localN; j++){
			printf("%f  ", f0[i][j]);
			//printf("%f  ", f1[i][j]);
		}
			printf("\n");
	}
	
	f0 = originalf0;
	f1 = originalf1;
	fend = originalfend;
	free(f0);
	free(f1);
	free(fend);
	
	MPI_Finalize();
	return 0;
}
	
void fInitOuterBounds(double **array, int localN){
	int i,j;
	for (int i=0; i<localN; i++){
		for(int j=0; j<localN; j++){
			array[0][j] = 0;
			array[localN-1][j] = 0;
			array[i][0] = 0;
			array[i][localN-1] = 0;
		//printf("%f  ", array[i][j]);
		}
		//printf("\n");
	}
}
