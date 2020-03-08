#include "demo4.h"
#include "array_util.h"

#include <mpi.h>

void main(int argc, char** argv)
{
    /* Data arrays */
    double *x;
    double s, total_sum;
    double mean;

    /* MPI variables */
    int my_rank, nprocs;
    int tag = 0;
    int count;
    int source;
    MPI_Status status;

    /* Other local variables */
    int i, p;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    set_rank(my_rank);  /* Needed to log printing */

    /* Total number of intervals */
    int num_intervals = pow2(12);


    if (my_rank == 0)
    {
        /* Compute integral */
        for(int p = 1; p < nprocs; p++)
        {
            int dest = p;
            /* Send interval information to each processor */
        }
    }
    else
    {
        int source = 0;
        /* Receive interval information from each processor */

    }

    /* Create array using interval */

    /* Aggregate sum over each interval on root node */

    MPI_Finalize();
}