#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>

void mpi_debug (void)
{
    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
    {
        printf("\n");
        printf("Getting setup for parallel debugging\n");
        printf("------------------------------------\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for(int i = 0; i < nprocs; i++)
    {
        if (rank == i)
        {
            printf("Proc %d with process %d is waiting to be attached\n",rank,getpid());
            fflush(stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    int ii = 0;
    while (ii == 0)  /* (gdb) set ii=1 */
    {
        /* Code will stop here;  set ii=1 to continue in gdb */
    }
}
