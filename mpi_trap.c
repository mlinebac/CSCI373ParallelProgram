// Author Matt Lineback
// CSCI 373 Homework03
// source code adapted from Program 3.1 in Pacheco
// 
#include <stdio.h>
#include <string.h>
#include <mpi.h>

const int MAX_STRING = 100;

double f(double x) {
  return x*x;
}

int main(int argz, char **argv) {
  
  FILE *fp;
  int comm_sz;
  int my_rank;
  
  MPI_Init(NULL,NULL);
  MPI_Comm_size(MPI_COMM_WORLD,&comm_sz);
  MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
  
  MPI_Barrier(MPI_COMM_WORLD);
  double local_start, local_finish, local_elapsed, elapsed;
  
  local_start = MPI_Wtime();
  double a=1.0;
  double b=3.0;
  int n;
  double h = (b-a)/(double)n;
  double sum = 0.0;
  if (my_rank == 0) {
    sum = (f(a)+f(b))/2.0;
  }

  for (int i=my_rank+1; i<=n-1; i=i+comm_sz) {
    double x_i = a + i*h;
    sum += f(x_i);
  }
  double solution = h*sum;

  if (my_rank !=0) {
    MPI_Send(&solution,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
  } else {
    for (int q=1; q<comm_sz; ++q) {
      double other_sum;
      MPI_Recv(&other_sum,1,MPI_DOUBLE,q,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      solution += other_sum;
    }

  fp = fopen("timing_data05.txt", "a");
  if (fp == NULL) {
    printf("Can't open this .txt file\n"); 
    return 1;
  } else {
     fprintf(fp,"%s %10s %10s %10s\n","prob_Size", "Num_proc", "Proc_time", "Solution");
     fprintf(fp,"%d %d %f %0.13f\n", n, comm_sz, (MPI_Wtime() - local_start), solution);
     fclose(fp);
  }
 }
  MPI_Finalize();
  return 0;
}

