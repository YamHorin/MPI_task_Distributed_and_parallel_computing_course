#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


//last change 1.8 13:56
    
enum tags {WORK, STOP ,DONE};
const int ROOT = 0;

double heron(int S, double epsilon, int *counter);
int test(int n_start, int n_end,int rank ,  double epsilon, int *min_number_max_counter);

int main(int argc, char **argv)
{
    int my_rank, num_procs , counter;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs); 

    double t_start = MPI_Wtime();

    int N = 100;
    double epsilon =0.01;
    int chunk_size = 10;
    switch (argc)
    {
    case 2:
        N = atoi(argv[1]);
        break;
    case 3:
         N = atoi(argv[1]);
        epsilon = atof(argv[2]);
        break;
    case 4:
         N = atoi(argv[1]);
        epsilon = atof(argv[2]);
        chunk_size = atoi(argv[3]);
        break;
    }
    if (chunk_size >N)
        chunk_size = N;
	if (num_procs < 2)
	{
		printf("eroor number of processes must be bigger then 2\n");
		MPI_Finalize();

	}
    int min_num = N;
    int data_pair[2];
    if (my_rank==0) //master
        {
            MPI_Status status;
            data_pair[1] =0;
            data_pair[0]  = 0;
            // each "task"  consists of finding max of CHUNK_SIZE numbers
            int tasks_sent = 0;
            // send  some work to each worker
            for(int worker_rank = 1; worker_rank < num_procs; worker_rank++){
                int start;
                start = (1 + ((worker_rank-1) * chunk_size));
                if (start<N)
                {
                    MPI_Send(&start,1, MPI_INT,worker_rank,WORK,MPI_COMM_WORLD );
                    tasks_sent++;
                }
                else
                    MPI_Send(&start,1, MPI_INT,worker_rank,STOP,MPI_COMM_WORLD );
            }
            // recv and send more work
            int tasks = N/chunk_size; // total number of "tasks". We assume remainder is 0
            if (N%chunk_size)
                tasks++;
            for(int tasks_done=0; tasks_done < tasks; tasks_done++)
            {
                int localMax;
                MPI_Recv(&localMax, 1, MPI_INT, MPI_ANY_SOURCE,
		                    DONE, MPI_COMM_WORLD, &status);      
                int tasks_not_sent_yet = tasks - tasks_sent;
                if (tasks_not_sent_yet > 0) {
                    int start;
                    start = 1 + (((tasks_sent-1) * chunk_size));
                    MPI_Send(&start,1, MPI_INT,status.MPI_SOURCE,WORK,MPI_COMM_WORLD );
                    tasks_sent++;
                }
                else {
                    /* send STOP message. message has no data */
                    int dummy;
                    MPI_Send(&dummy,0, MPI_INT, status.MPI_SOURCE,
                            STOP, MPI_COMM_WORLD);
                }
            }
        }
    else //worker
    {
        MPI_Status status;
        int tag;
        do
        {
            int j , finish;
            MPI_Recv(&j , 1 , MPI_INT , ROOT,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
            finish =j+chunk_size; 
             tag = status.MPI_TAG;
            if (finish>N)
                finish = N;
            if (tag==WORK)
            {
                counter  = test(j, finish ,my_rank , epsilon ,&min_num);
                if (data_pair[0]<counter){
                    data_pair[1] = min_num;//return heron func
                    data_pair[0] = counter;//counter
            }
            MPI_Send(&counter,1,MPI_INT,ROOT,
			         DONE,MPI_COMM_WORLD);
            }
            
        } while (tag != STOP);  
     }
    int reduction_result[2];
    MPI_Reduce(data_pair, reduction_result, 1, MPI_2INT, MPI_MAXLOC, 0, MPI_COMM_WORLD);
    if(my_rank == 0)
    {
        printf("\nnumber requiring max number of iterations : %d (number of iterations : %d)\n",reduction_result[1],reduction_result[0]);
        printf("sequential time: %f secs\n", MPI_Wtime() - t_start);    
    }
    MPI_Finalize();
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
	
    //printf("rank = %d , max_counter = %d , min_num = %d\n\n",rank , max_counter , *min_number_max_counter);
	return max_counter;
}
