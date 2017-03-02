#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>
//#ifndef M_PI
//#define M_PI 3.14159265359

void writeheader(int N, int end) {
	FILE *fp;
	fp = fopen("testoutfile.pgm", "w");
	if (fp == NULL) {
		printf("sorry can't open testoutfile.pgm. Terminating.\n");
		exit(1);
	}
	else {
		// print a table header
		fprintf(fp, "%s\n%d %d\n%s\n", "P2", N, end, "255");
		fclose(fp);
	}
}

void writerow(int N, double **rawdata) {
	FILE *fp;
	fp = fopen("testoutfile.pgm", "a");
	if (fp == NULL) {
		printf("sorry can't open outfile.pgm. Terminating.\n");
		exit(1);
	}
	else {
		for (int i=0; i<N; i++) {
				for(int j=0; i<N-1; j++){
			int val = rawdata[i][j]*127+127;
			fprintf(fp,"%d ",val);
		}
		fprintf(fp,"\n");
	}
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

void fInitLeftOuterBounds(double **array, int localN);
void fInitRightOuterBounds(double **array, int localN);
void printArray(double **array, int localN);

int main(int argc, char *argv[]) {
	int comm_sz;
	int my_rank;
	int N = 6;
	int localN;
	int end = 6;//end=20N is roughly 1 period
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
	double **temp;
	
	double forOutput[N];
	
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
	
	for (int i=0; i<localN; i++){
		for (int j=0; j<=localN; j++){
			x = localx + (double)i*1.0/(N-1);
			y = localy + (double)j*1.0/(N-1);
			f0[i][j] = initialCondition(x,y);
			f1[i][j] = initialCondition(x,y);
	}
	
}	
//printArray(f1, localN);
	if (my_rank == 0){
		fInitLeftOuterBounds(f0,N);
		fInitLeftOuterBounds(f1,N);
		fInitLeftOuterBounds(fend,N);
	
	}
	
	if(my_rank == comm_sz-1){
		fInitRightOuterBounds(f0,localN);
		fInitRightOuterBounds(f1,localN);
		fInitRightOuterBounds(fend,localN);
	}
	//printArray(f1, localN);
	if (writeoutput){
		MPI_Gather(f1, localN, MPI_DOUBLE, forOutput, localN, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		if (my_rank == 0){
			writeheader(N, end);
			writerow(N, forOutput);
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
		MPI_Sendrecv(&f1[1][localN-1],1,MPI_DOUBLE,partnerRight,0,&leftneighbor,1,MPI_DOUBLE,partnerLeft,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		//send left
		MPI_Sendrecv(&f1[1][0],1,MPI_DOUBLE,partnerLeft,0,&rightneighbor,1,MPI_DOUBLE,partnerRight,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

		//put next step in f2:
		//compute interior of my domain
		for (int i=1; i<localN-1; i++){
		for (int j=1; j<=localN-2; j++){
			fend[i][j]= 0.01*(leftneighbor+f1[i+1][j]+f1[i][j-1]+f1[i][j+1]-4*f1[i][j])+2*f1[i][j]-f0[i][j];
			
		}
	}
		//printArray(fend, localN);
		//compute left edge of my domain
		if (my_rank!=0) {
			int i=0;
			fend[i][j] = 0.01*(f1[i-1][j]+rightneighbor+f1[i][j-1]+f1[i][j+1]-4*f1[i][j])+2*f1[i][j]-f0[i][j];
			
		}
		//compute right edge of my domain
		if (my_rank!=comm_sz-1) {
			int i=localN-1;
			fend[i][j] = 0.01*(f1[i-1][j]+f1[i+1][j]+f1[i][j-1]+f1[i][j+1]-4*f1[i][j])+2*f1[i][j]-f0[i][j];
			
		}
		printArray(fend, localN);
		//write output
		if (writeoutput) {
			MPI_Gather(fend,localN,MPI_DOUBLE,forOutput,localN,MPI_DOUBLE,0,MPI_COMM_WORLD);
			if (my_rank==0) {
				writerow(N,forOutput);
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
	free(f0);
	free(f1);
	free(fend);
	
	MPI_Finalize();
	return 0;
}
	
void fInitLeftOuterBounds(double **array, int localN){
	int i,j;
	for (int i=0; i<localN; i++){
		for(int j=0; j<localN; j++){
			array[0][j] = 0;
			array[i][0] = 0;
		}
	}
}
void fInitRightOuterBounds(double **array, int localN){
	int i,j;
	for (int i=0; i<localN; i++){
		for(int j=0; j<localN; j++){
			array[localN-1][j] = 0;
			array[i][localN-1] = 0;
		}
	}
}
void printArray(double **array, int localN){
	int i,j;
	for (int i=0; i<localN; i++){
		for(int j=0; j<localN; j++){
			printf("%f  ", array[i][j]);
		}
			printf("\n");
	}
}
