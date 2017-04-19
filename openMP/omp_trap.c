#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

double f(double);
void TrapParallel(double a, double b, int n, double* global_result_p);

int main(int argc, char* argv[]){
	
	int n;
	double a = 1;
	double b = 3;
	int thread_count;
	double global_result = 0.0;
	double start;
	double end;
	double timeDiff;
	
	thread_count = strtol(argv[1], NULL, 10);
	printf("Enter number of trapezoids\n");
	scanf("%d", &n);
	
	if (n % thread_count !=0) {
		fprintf(stderr, "n must be evenly divisible by thread_count/n");
		exit(0);
	}
	
	start = omp_get_wtime();
#	pragma omp parallel num_threads(thread_count)
	TrapParallel(a, b, n, &global_result);
	end = omp_get_wtime();
	printf("With n = %d trapezoids, our estimate\n", n);
	printf("of the integral from %f to %f = %.14e\n", a, b, global_result);
	timeDiff = (end-start);
	printf("Parallel time is %f\n", timeDiff);
	
	return 0;
	
}
	
double f(double x) {
		return x*x;
	}
	
void TrapParallel(double a, double b, int n, double* global_result_p){
	double h, x, my_result;
	double local_a, local_b;
	int i, local_n;
	int my_rank = omp_get_thread_num();
	int thread_count = omp_get_num_threads();
	
	h = (b-a)/(double)n;
	local_n = n/thread_count;
	local_a = a + my_rank*local_n*h;
	local_b = local_a + local_n*h;
	my_result = (f(local_a) + f(local_b))/2.0;
	
	for(i = 1; i <= local_n-1; i++){
		x = local_a + i*h;
		my_result += f(x);
	}
	my_result = my_result*h;
	//printf("result from thread %d = %f\n", my_rank, my_result);
	
#	pragma omp critical 
	*global_result_p += my_result;
}
	
