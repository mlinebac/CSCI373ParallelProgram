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
void printArray(double *array, int length) {
  printf("[");
  int i;
  for (i = 0; i < length; i++) {
		if (i != 0) {
      printf("  ");
    }
    printf("%g", array[i]);
  }
  printf("]");
}

	
int main(int argc, char *argv[]) {
	int comm_sz;
	int my_rank;
	int N = 6;
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
	
	//double * f = (double *)malloc(localN*sizeof(double));
	//double * f1 = (double *)malloc(localN*sizeof(double));
	//double * f2 = (double *)malloc(localN*sizeof(double));
	double localx = 1.0/(N-1)*my_rank*localN;
	double localy = 1.0/(N-1)*my_rank*localN;
	double f[6][6];
	double x;
	double y;
	
	for (int i=0; i<localN; i++){
		x = localx + (double)i*1.0/(N-1);
		y = localy + (double)i*1.0/(N-1);
		printf("this is x%f\t",x);
		printf("this is y%f\t\n",y);
}
	int j;
	for (int i=0; i<localN; i++){
		for (int j=0; j<=localN-1; j++){
		f[i][j] = initialCondition(x,y);
		printf("%f\t", f[i][j]);
		/*f[0][j] = 0;
		f[1][j] = 0;
		f[i][0] = 0;
		f[i][1] = 0;
		*/
	}
	printf("\n");
}	
	
	//f[i][j] = initialCondition(x,y);
	MPI_Finalize();
	return 0;
}
	
