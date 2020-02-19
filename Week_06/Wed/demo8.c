#include <stdio.h>
#include <mpi.h>

#include "array_util.h"

void main(int argc, char** argv)
{
    MPI_Request request1;
    MPI_Request request2;

    MPI_Init(&argc, &argv);

    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    int tag = 0;
    random_seed();

    int root = 0;
    double x = random_number() + my_rank;                
    printf("My x is %f\n",x);
    MPI_Bcast(&x,1,MPI_DOUBLE,root,MPI_COMM_WORLD);                
    printf("Now my x is %f\n",x);

    MPI_Finalize();
}