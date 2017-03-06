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
			//for(int i=0; i<end; i++){
			for(int j=0; j<N; j++){
				int val = rawdata[i][j]*127;
				fprintf(fp,"%d ", val);
			}
			fprintf(fp,"\n");
			
		//}
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

void fInitLeftOuterBounds(double **array, int localN, int N);
void fInitRightOuterBounds(double **array, int localN, int N);
void printArray(double **array, int localN, int N);

int main(int argc, char *argv[]) {
	int comm_sz;
	int my_rank;
	int N = 12;
	int localN;
	int end = 12;//end=20N is roughly 1 period
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
	double **forOutput = malloc(N * sizeof *forOutput);
	double **temp;
	
	int i,j;
	for(i=0; i<N; i++){
		 f0[i] = malloc(N*sizeof *f0[i]);
		 f1[i] = malloc(N*sizeof *f1[i]);
		 fend[i] = malloc(N*sizeof *fend[i]);
		 forOutput[i] = malloc(N*sizeof *forOutput[i]);	 
		 
	}
	
	double localx = 1.0/(N-1)*my_rank*localN;
	double localy = 1.0/(N-1)*my_rank*localN;
	
	double **originalf0 = f0;
	double **originalf1 = f1;
	double **originalfend = fend;
	double **originalforOutput = forOutput;
	
	
	double x;
	double y; 
	
	for (int i=1; i<end-1; i++){
		for (int j=1; j<N-1; j++){
			x = localx + (double)i*1.0/(N-1);
			y = localy + (double)j*1.0/(N-1);
			f0[i][j] = initialCondition(x,y);
			f1[i][j] = initialCondition(x,y);
			
		}
	}	
	
	/*if (my_rank == 0){
		fInitLeftOuterBounds(f0,localN, N);
		fInitLeftOuterBounds(f1,localN, N);
		fInitLeftOuterBounds(fend, localN, N);
	
	}
	
	if(my_rank == comm_sz-1){
		fInitRightOuterBounds(f0,localN, N);
		fInitRightOuterBounds(f1,localN, N);
		fInitRightOuterBounds(fend, localN, N);
		
	}
*/
	//printArray(f1, localN, N);
	if (writeoutput){
		MPI_Gather(f1, localN, MPI_DOUBLE, forOutput, localN, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		if (my_rank == 0){
			writeheader(localN, N);
			writerow(localN, N, forOutput);
			//printArray(forOutput,localN, N);
		}
	}
		
	double leftneighbor=0.0;
	double rightneighbor=0.0;
	
	int partnerRight = my_rank+1;
	int partnerLeft = my_rank-1;
	if (partnerRight==comm_sz) {
		partnerRight=MPI_PROC_NULL;
	}
	if (partnerLeft==-1) {
		partnerLeft=MPI_PROC_NULL;
	}

	//main loop
	int step = 2;//current step
	while (step<=end) {
		//send right
		MPI_Sendrecv(&f1[0][N],1,MPI_DOUBLE,partnerRight,0,&leftneighbor,1,MPI_DOUBLE,partnerLeft,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		//send left
		MPI_Sendrecv(&f1[localN-1][N],1,MPI_DOUBLE,partnerLeft,0,&rightneighbor,1,MPI_DOUBLE,partnerRight,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

		//put next step in f2:
		//compute interior of my domain
		for (int i=1; i<localN-1; i++){
			for (int j=1; j<N-1; j++){
				fend[i][j]= 0.01*(f1[i-1][j]+f1[i+1][j]+f1[i][j-1]+f1[i][j+1]-4*f1[i][j])+2*f1[i][j]-f0[i][j];
			
			
			}
		}
		//printArray(fend, localN, N);
		
		//compute left edge of my domain
		if (my_rank!=0) {
			for (int i=1; i<localN-1; i++){
				for (int j=1; j<N-1; j++){
					fend[i][j] = 0.01*(f1[i-1][j]+leftneighbor+f1[i][j-1]+f1[i][j+1]-4*f1[i][j])+2*f1[i][j]-f0[i][j];
				}
			}
		}
		//compute right edge of my domain
		if (my_rank!=comm_sz-1) {
			for (int i=1; i<localN-1; i++){
				for (int j=1; j<N-1; j++){
					fend[i][j] = 0.01*(rightneighbor+f1[i+1][j]+f1[i][j-1]+f1[i][j+1]-4*f1[i][j])+2*f1[i][j]-f0[i][j];
				}	
			}
		}
		//printArray(fend, localN);
		//write output
		if (writeoutput) {
			MPI_Gather(fend,localN,MPI_DOUBLE,forOutput,localN,MPI_DOUBLE,0,MPI_COMM_WORLD);
			if (my_rank==0) {
				
				writerow(N, localN,forOutput);
				printArray(forOutput, localN, N);
			}
		}

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
	
void fInitLeftOuterBounds(double **array, int localN, int N){
	int i,j;
	for (int i=0; i<localN; i++){
		for(int j=0; j<N; j++){
			array[0][j] = 0;
			array[localN-1][j];
			array[i][0] = 0;
			array[0][N-1]=0;
		}
	}
}
void fInitRightOuterBounds(double **array, int localN, int N){
	int i,j;
	for (int i=0; i<localN; i++){
		for(int j=0; j<N; j++){
			array[localN-1][j] = 0;
			array[i][N-1] = 0;
			array[i][0] = 0;
			array[0][j] = 0;
		}
	}
}
void printArray(double **array, int localN, int N){
	int i,j;
	for (int i=0; i<localN; i++){
		for(int j=0; j<N; j++){
			printf("%f  ", array[i][j]);
		}
			printf("\n");
	}
}

