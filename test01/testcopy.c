#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>
//#ifndef M_PI
//#define M_PI 3.14159265359

void writeheader(int N, int end) {
	FILE *fp;
	char outputfile[20];
	int k;
	for (int k=0; k<end; k++){
		sprintf(outputfile,"output%04d.pgm", k);
		fp = fopen(outputfile, "w");
		
		if (fp == NULL) {
		printf("sorry can't open output.pgm. Terminating.\n");
		exit(1);
		}
		else {
		fprintf(fp, "%s\n%d %d\n%s\n", "P2", N, end, "255");
		fclose(fp);
		}
	}
}
void writerow(int end, int N, double **rawdata) {
	FILE *fp;
	char outputfile[20];
	int i,j;
	for(int i=0; i<end; i++){
		sprintf(outputfile,"output%04d.pgm",i);
		fp = fopen(outputfile,"a");
	
		if (fp == NULL) {
			printf("sorry can't open outfile.pgm. Terminating.\n");
			exit(1);
		}
		else {
			for(int j=0; j<N; j++){
				int val = rawdata[i][j]*127;
				fprintf(fp,"%d ", val);
			}
			fprintf(fp,"\n");
			fclose(fp);	
		}
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

void fInitLeftOuterBounds(double **array, int localEnd, int N);
void fInitRightOuterBounds(double **array, int localEnd, int N);
void printArray(double **array, int localEnd, int N);

int main(int argc, char *argv[]) {
	int comm_sz;
	int my_rank;
	int N = 8;
	int localEnd;
	int end = 8;//end=20N is roughly 1 period
	int writeoutput = 1;//0 for false
	int i,j;
	
	MPI_Init(NULL,NULL);
	MPI_Comm_size(MPI_COMM_WORLD,&comm_sz);
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
	
	localEnd = N/comm_sz;//assuming divisible
	
	if (N%comm_sz!=0) {
		printf("ERROR N must be divisible by comm_sz. Terminating.");
		exit(1);
	}
	if (localEnd<2) {
		printf("ERROR N (%d) must be at least twice as large as number of processes (%d). Terminating.",N,comm_sz);
		exit(1);
	}
	
	double **f0 = malloc(N * sizeof(double));
	double **f1 = malloc(N * sizeof (double));
	double **fend = malloc(N * sizeof (double));
	double **forOutput = malloc(N * sizeof (double));
	double **temp;
	
	for(i=0; i<N; i++){
		 f0[i] = malloc(N * sizeof (double));
		 f1[i] = malloc(N * sizeof (double));
		 fend[i] = malloc(N * sizeof (double));
		 forOutput[i] = malloc(N * sizeof (double));	 
		 
	}
	
	double localx = 1.0/(N-1)*my_rank*localEnd;
	double localy = 1.0/(N-1)*my_rank*localEnd;
	
	double **originalf0 = f0;
	double **originalf1 = f1;
	double **originalfend = fend;
	double **originalforOutput = forOutput;
	
	
	double x;
	double y; 
	
	for (int i=1; i<N-1; i++){
		for (int j=1; j<N-1; j++){
			x = localx + (double)i*1.0/(N-1);
			y = localy + (double)j*1.0/(N-1);
			f0[i][j] = initialCondition(x,y);
			f1[i][j] = initialCondition(x,y);
			
		}
	}	
	
	/*printArray(f1, localEnd, N);
	for(int j=0; j<N; j++){
		printf("%f ",f1[localEnd/2-1][j]);
		
	}
	printf("\n");
	*/
	if (writeoutput){
		MPI_Gather(f1, N, MPI_DOUBLE, forOutput, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		if (my_rank == 0){
			writeheader(localEnd, N);
			writerow(localEnd, N, forOutput);
			//printf("from proc 0 \n");
			printArray(forOutput,end, N);
			
		}
	}
	
	double upNeighbor[N];
	double downNeighbor[N];
	
	 
	int partnerDown = my_rank+1; 
	int partnerUp = my_rank-1; 
	
	if (partnerDown == comm_sz) {
		partnerDown=MPI_PROC_NULL;
	}
	if (partnerUp == -1) {
		partnerUp=MPI_PROC_NULL;
	}
	
	//main loop
	int step = 2;//current step
	while (step<=end) {
		//send bottom
		MPI_Sendrecv(&f1[localEnd-1],N,MPI_DOUBLE,partnerDown,0,&downNeighbor,N,MPI_DOUBLE,partnerUp,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		//send top
		MPI_Sendrecv(&f1[localEnd],N,MPI_DOUBLE,partnerUp,0,&upNeighbor,N,MPI_DOUBLE,partnerDown,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		
		
		//put next step in f2:
		//compute interior of my domain
		if (my_rank!=0){
			for (int i=localEnd; i<N-2; i++){
				for (int j=1; j<N-2; j++){
					fend[i][j]= 0.01*(f1[i-1][j]+f1[i+1][j]+f1[i][j-1]+f1[i][j+1]-4*f1[i][j])+2*f1[i][j]-f0[i][j];
				}
			}
		}
		else{
			for (int i=1; i<localEnd-1; i++){
				for (int j=1; j<N-2; j++){
					fend[i][j]= 0.01*(f1[i-1][j]+f1[i+1][j]+f1[i][j-1]+f1[i][j+1]-4*f1[i][j])+2*f1[i][j]-f0[i][j];
				}
			}
		}
		//printArray(fend, localEnd, N);
		//printf("\n");
		//compute top edge of my domain
		if (my_rank!=0) {
			//int i = localEnd; 
			for (int j=1; j<N-2; j++){
					fend[4][j] = 0.01*(upNeighbor[j]+f1[i+1][j]+f1[i][j-1]+f1[i][j+1]-4*f1[i][j])+2*f1[i][j]-f0[i][j];
			}
		}
		
		//printArray(fend, localN, N);
		//compute bottom edge of my domain
		if (my_rank!=comm_sz-1) {
			//int i = localEnd-1;
			for (int j=1; j<N-2; j++){
					fend[3][j] = 0.01*(f1[i-1][j]+downNeighbor[j]+f1[i][j-1]+f1[i][j+1]-4*f1[i][j])+2*f1[i][j]-f0[i][j];
			}	
		}
		//printArray(fend, localN, N);
		//write output
		if (writeoutput) {
			MPI_Gather(fend,N,MPI_DOUBLE,forOutput,N,MPI_DOUBLE,0,MPI_COMM_WORLD);
			if (my_rank==0) {
				writerow(N, localEnd,forOutput);
				//printf("from proc 0 after end steps\n");
				//printArray(forOutput, end, N);
			}
		}
		//printArray(forOutput, localN, N);
		//rearrange pointers for next step
		temp=f0;
		f0=f1;
		f1=fend;
		fend=temp;

		step++;

	}

	f0 = originalf0;
	f1 = originalf1;
	fend = originalfend;
	forOutput = originalforOutput;
	
	free(f0);
	free(f1);
	free(fend);
	free(forOutput);
	
	
	MPI_Finalize();
	return 0;
}
	
void fInitLeftOuterBounds(double **array, int localEnd, int N){
	int i,j;
	for (int i=0; i<localEnd; i++){
		for(int j=0; j<N; j++){
			array[0][j] = 0;
			array[localEnd-1][j];
			array[i][0] = 0;
			array[0][N-1]=0;
		}
	}
}
void fInitRightOuterBounds(double **array, int localEnd, int N){
	int i,j;
	for (int i=0; i<localEnd; i++){
		for(int j=0; j<N; j++){
			array[localEnd-1][j] = 0;
			array[i][N-1] = 0;
			array[i][0] = 0;
			array[0][j] = 0;
		}
	}
}
void printArray(double **array, int localEnd, int N){
	int i,j;
	for (int i=0; i<localEnd; i++){
		for(int j=0; j<N; j++){
			printf("%f  ", array[i][j]);
		}
			printf("\n");
	}
}


