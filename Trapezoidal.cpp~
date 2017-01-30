// Trapezoidal.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

	double f(double);
	double area(double,double,int,double);

	int main(void) {

		double n = 1000, local_n;
		double a = 1,b = 3, local_a, local_b;
		double h;
		double local_int, total_int;
		int source;

		int comm_sz;
		int my_rank;

		MPI_Init(NULL, NULL);
		MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
		MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

		h = (b - a) / n;//this is the delta x or the step size of our approximation
		local_n = n / comm_sz;

		local_a = a + my_rank * local_n * h;
		local_b = local_a + local_n * h;
		local_int = area(local_a, local_b, local_n, h);

		if (my_rank != 0) {
			MPI_Send(&local_int, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		}
		else {
			total_int = local_int;
			for (source = 1; source < comm_sz; source++) {
				MPI_Recv(&local_int, 1, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				total_int += local_int;
			}
		}
		if (my_rank == 0) {
			printf("With n = %d trapezoids, our approximation\n", n);
			printf("of the integral from %f to %f = %.15e\n", a, b, total_int);
		}
		MPI_Finalize();
		return 0;
	}

	double area(double a, double b, int trap_count, double h) {

		double approx, tempX;
		int i;

		approx = (f(a) + f(b)) / 2.0;
		for (i = 1; i < trap_count-1; i++) {
			tempX = a + i*h;
			approx += f(tempX);
		}
		approx = approx * h;
		return approx;
	}

	double f(double x) {
		return x*x;
	}
