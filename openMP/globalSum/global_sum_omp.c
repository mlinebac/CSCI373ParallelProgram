/*
 * Author: Matt Lineback
 * Global Sum Example for OpenMP
 * 4/29/2017
*/


#include <stdio.h>
#include <stdlib.h>
#include <omp.h>


int main(int argc, char* argv[]) {
	int n, i, thread_count;
	thread_count = strtol(argv[1], NULL, 10);
	n = strtol(argv[2], NULL, 10);
    double global_sum = 0.0;
    double start, end, timediff;
	start = omp_get_wtime();
    
    # pragma omp parallel for num_threads(thread_count) reduction(+: global_sum)
       for(i=0; i<n; i++){
	      global_sum += i;
	   }
    
    end = omp_get_wtime();
	timediff = end-start;
	printf("time elapsed %f\n", timediff);
    printf("global sum = %f\n", global_sum);
    
    return 0;
}
