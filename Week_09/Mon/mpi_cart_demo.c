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

    double qmem[3];    
    double *q = &qmem[1];

    q[-1] = -999.;
    q[0] = rank;
    q[1] = 999;

    /* Fill in q[-1] and q[1] with data from neighbors */
    MPI_Comm comm_cart;
    int ndim = 1;
    int dims[1] = {nprocs};
    int periodicity[1] = {1};
    int reorder = 0;
    MPI_Cart_create(MPI_COMM_WORLD, ndim, dims,  periodicity, reorder, &comm_cart);

    /* Exchange left */
    int dir = 0;
    int disp = 1;  
    int source, dest;
    int ierr = MPI_Cart_shift(comm_cart, dir, disp, &source, &dest);

    int tag = 0;
    MPI_Sendrecv(&q[0], 1, MPI_DOUBLE, dest, tag, 
                            &q[-1], 1, MPI_DOUBLE, source, tag, comm_cart, MPI_STATUS_IGNORE);

    disp = -1;  /* exchange right */
    MPI_Cart_shift(comm_cart, dir, disp, &source, &dest);

    tag = 1;
    ierr = MPI_Sendrecv(&q[0], 1, MPI_DOUBLE, dest, tag, 
                        &q[1], 1, MPI_DOUBLE, source, tag, comm_cart, MPI_STATUS_IGNORE);

    printf("Rank %d : (%g, %g, %g)\n",rank, q[-1], q[0], q[1]);

    MPI_Finalize();
}

