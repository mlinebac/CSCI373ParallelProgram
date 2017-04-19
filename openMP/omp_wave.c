#include <stdio.h>
#include <string.h>
#include <omp.h>
#include <stdlib.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


int main(int argc, char *argv[]) {
	int N = 6;
	int localN;
	//int end = 6;//end=20N is roughly 1 period
	//int thread_count = 4;
	//int my_rank = omp_get_thread_num();
	//int comm_sz = omp_get_num_threads();
	//localN = N/my_rank;
	
	double * f0 = (double *)malloc(localN*sizeof(double));
	double * f1 = (double *)malloc(localN*sizeof(double));
	double * f2 = (double *)malloc(localN*sizeof(double));
	//double * temp;// = (double *)malloc(localN*sizeof(double));
	
	//keep original pointers so we can free later
	double * originalf0 = f0;
	double * originalf1 = f1;
	double * originalf2 = f2;
	
	double localx0 = 1.0/(N-1)*my_rank*localN;
	
	for (int i=0; i<localN; i++) {
		double x = localx0 + (double)i*1.0/(N-1);
		f0[i] = sin(M_PI*x);
		f1[i] = sin(M_PI*x);
	}
	if (my_rank==0) {
		f0[0] = 0;
		f1[0] = 0;
		f2[0] = 0;
	}
	if (my_rank==comm_sz-1) {
		f0[localN-1] = 0;
		f1[localN-1] = 0;
		f2[localN-1] = 0;
	}
	
	for (int i=1; i<localN-1; i++) {
			f2[i] = 0.01*(f1[i-1]-2*f1[i]+f1[i+1]) + 2*f1[i] - f0[i];
		}
	/*double leftneighbor = f1[i-1];
	double rightneighbor = f1[i+1];
	int step = 2;
	while (step<=end){
		
# 	pragma omp parallel for num_threads(thread_count) private(leftneighbor, rightneighbor)
		
		for (int i=1; i<localN-1; i++) {
			f2[i] = 0.01*(f1[i-1]-2*f1[i]+f1[i+1]) + 2*f1[i] - f0[i];
		}
		for(int i = 0; i<N; i++){
	printf("%f\n", f2[i]);
}
		//compute left edge of my domain
		if (my_rank!=0) {
			int i=0;
			f2[i] = 0.01*(leftneighbor-2*f1[i]+f1[i+1]) + 2*f1[i] - f0[i];
		}

		//compute right edge of my domain
		if (my_rank!=comm_sz-1) {
			int i=localN-1;
			f2[i] = 0.01*(f1[i-1]-2*f1[i]+rightneighbor) + 2*f1[i] - f0[i];
		}

		//rearrange pointers for next step
		temp=f0;
		f0=f1;
		f1=f2;
		f2=temp;

		step++;
	}
	*/
	f0 = originalf0;
	f1 = originalf1;
	f2 = originalf2;
	free(f0);
	free(f1);
	free(f2);



	return 0;
}
