//Author Matt Lineback
//CSCI 373 Parallel Programming 
//HW05

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>


void print_file_header();
void print_file_data();
void print_data(long int number_of_tosses, int thread_count, double timediff, double pi_estimate);

int main(int argc, char* argv[]){
	
	long int number_of_tosses, number_in_circle;
	int toss, thread_count;
	double x, y, distance_squared, start, end, timediff;
	unsigned short xi[3];
	
	if (argv[1] == NULL){
		printf("enter number of threads as first argument\n");
		exit(1);
	}else{
		thread_count = strtol(argv[1], NULL, 10);
		}
	
	if (argv[2] == NULL){
		printf("enter number of tosses as second argument\n");
		exit(1);
	}else{
		number_of_tosses = strtoll(argv[2], NULL, 10);
	}
	start = omp_get_wtime();
	
	xi[0]=1;
	xi[1]=1;
	xi[2]=omp_get_thread_num();
	number_in_circle = 0;
	
#	pragma omp parallel for num_threads(thread_count) \
		reduction(+:number_in_circle) \
		firstprivate(xi) private(x, y, distance_squared) 
	for (toss = 0; toss < number_of_tosses; toss++) {
		x = erand48(xi);
		y = erand48(xi);
		distance_squared = x*x + y*y;
		if (distance_squared <= 1){
			number_in_circle++;
		}
	}
	
	double pi_estimate = 4*number_in_circle/((double) number_of_tosses);
	printf("estimate of pi %.14f\n", pi_estimate);
	end = omp_get_wtime();
	timediff = end-start;
	printf("time elapsed %f\n", timediff);
	//print_file_header();
	print_file_data(number_of_tosses, thread_count, timediff, pi_estimate);

	return 0;
}

void print_file_header(){
	FILE *fp;
		fp = fopen("monte_carlo_data.txt", "w");
		if (fp == NULL) {
			printf("sorry can't open output file\n");
			exit(1);
		}
		else {
			fprintf(fp,"%s %s %s %s\n","Tosses", "Threads", "TimeElapsed", "PiEstimation");
			fclose(fp);
		}
	}
	
void print_file_data(long int number_of_tosses, int thread_count, double timediff, double pi_estimate){
	FILE *fp;
		fp = fopen("monte_carlo_data.txt", "a");
		if (fp == NULL) {
			printf("sorry can't open output file\n");
			exit(1);
		}
		else {
			fprintf(fp,"%ld %d %f %0.14f\n", number_of_tosses, thread_count, timediff, pi_estimate);
			fclose(fp);
		}
	}
	
	
