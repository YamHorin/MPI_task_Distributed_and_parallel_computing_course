#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>



//last change 31.7 22:23
double heron(int S, double epsilon, int *counter);
int test(int n_start, int n_end,int rank ,  double epsilon, int *min_number_max_counter);

int main(int argc, char **argv)
{
	int my_rank, num_procs, counter;

	MPI_Init(&argc, &argv);
	double t_start = MPI_Wtime();

	int N = 100;
	double epsilon = 0.01;
	if (argc >= 2)
		N = atoi(argv[1]);
	if (argc == 3)
		epsilon = atof(argv[2]);

	int min_num = N;
	int data_pair[2];
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	if (my_rank==0)
	{
		if(N%(num_procs - 1)) //check if the num of process does npt divid
		{
		printf("the number of processes don't dividde by the number of processes (not include the master), exit program\n");
		MPI_Finalize();
		return 1;
		}
	}
	// Compute the workload division among non-master processes (ranks != 0)
	
    int n_per_process = (N + num_procs - 2) / (num_procs - 1);
    int n_start = 1 + (my_rank - 1) * n_per_process;
    int n_end = my_rank * n_per_process;
	counter = test(n_start, n_end ,my_rank, epsilon, &min_num);

	data_pair[1] = min_num; // return heron func
	data_pair[0] = counter; // counter
	int reduction_result[2];
	MPI_Reduce(data_pair, reduction_result, 1, MPI_2INT, MPI_MAXLOC, 0, MPI_COMM_WORLD);
	if (my_rank == 0)
	{
		printf("\nnumber requiring max number of iterations : %d (number of iterations : %d)\n", reduction_result[1], reduction_result[0]);
		printf("sequential time: %f secs\n", MPI_Wtime() - t_start);
	}
	MPI_Finalize();
	return 0;
}

double heron(int S, double epsilon, int *counter)
{
	double x = 10.0; // first approximation (we could do better ...)

	*counter = 0;
	double next_x = x;
	do
	{
		x = next_x;
		next_x = (x + S / x) / 2;
		(*counter)++; //  count number of iterations
	} while (fabs(next_x - x) > epsilon);

	return next_x;
}

int test(int n_start, int n_end,int rank ,  double epsilon, int *min_number_max_counter)
{
	int max_counter = 0;
	int count = 0;
	int result;
	if (rank!=0)
	{
		for (int i = n_start; i <= n_end; i++)
		{
		result = heron(i, epsilon, &count);
		if (count > max_counter)
			{
				max_counter = count;
				*min_number_max_counter = i;
			}
		count = 0;
		}
	}
	return max_counter;
}
