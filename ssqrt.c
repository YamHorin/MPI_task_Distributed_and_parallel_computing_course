#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

//last change 31.7 22:23
int test(int n_start, int n_end,  double epsilon, int *min_number_max_counter);
double heron(int S, double epsilon, int *counter);

int main(int argc, char **argv)
{
    int num_procs;
    MPI_Init(&argc, &argv); //just for the time function
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs); 
    if (num_procs!=1)
    {
        printf("the normal version can only work with 1 process , exit\n");
        MPI_Finalize();
        return 1;

    }
    int N = (argc > 1) ? atoi(argv[1]) : 100;
    double epsilon = (argc > 2) ? atof(argv[2]) : 0.01;
    double t_start = MPI_Wtime();
    int min_num;
	int counter = test(1, N , epsilon, &min_num);	
    printf("\nnumber requiring max number of iterations : %d (number of iterations : %d)\n", min_num, counter);
	printf("sequential time: %f secs\n", MPI_Wtime() - t_start);
	MPI_Finalize();
    return 0;

}

double heron(int S, double epsilon, int *counter)
{
    double x = 10.0; // first approximation (we could do better ...)

    *counter = 0;
	double next_x = x; 
    do {
		x = next_x;
        next_x = (x + S/x)/2;
		(*counter)++;   //  count number of iterations
    } while (fabs(next_x - x) > epsilon);
 
    return next_x; 
}
int test(int n_start, int n_end,  double epsilon, int *min_number_max_counter)
{
	int max_counter = 0;
	int count = 0;
	int result;
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
	return max_counter;
}