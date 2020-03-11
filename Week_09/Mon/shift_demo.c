#include <mpi.h>
#include <math.h>

#include <stdio.h>
#include <stdlib.h>

void main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    int myval = rank;

    MPI_Comm comm_cart;
    int ndim = 1;
    int dims[1] = {nprocs};
    int periodicity[1] = {0};
    int reorder = 0;
    MPI_Cart_create(MPI_COMM_WORLD, ndim, dims,  periodicity, reorder, &comm_cart);

    int dir = 0;
    int disp = 1;  /* Get a value from right neighbor */
    int source, dest;
    /* source : source of data for this proc
       dest   : Where this proc will send data 
    */
    MPI_Cart_shift(comm_cart, dir, disp, &source, &dest);
    if (source == MPI_PROC_NULL && disp == 1)
    {
        printf("Rank %d : No left neighbor found\n",rank);
    }
    else if (source == MPI_PROC_NULL && disp == -1)
    {
        printf("Rank %d : No right neighbor found\n",rank);
    }

    int tag = 0;
    MPI_Sendrecv_replace(&myval, 1, MPI_INTEGER, dest, tag, 
                                    source, tag, comm_cart, MPI_STATUS_IGNORE);

    printf("Rank %d : myval = %d (%d, %d)\n",rank, myval, source, dest);

    MPI_Finalize();
}

