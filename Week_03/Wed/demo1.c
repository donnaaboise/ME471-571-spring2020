#include <stdio.h>
#include <mpi.h>
#include "array_util.h"

void main(int argc, char** argv)
{
    int my_rank;
    int p;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    MPI_Comm_size(MPI_COMM_WORLD, &p);

    printf("My rank is %d\n", my_rank);

    /* 
    sleep(5.0);
    
    printf("Rank %d is done!\n",my_rank);    
    */

    MPI_Finalize();

}