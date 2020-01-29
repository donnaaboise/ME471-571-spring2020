#include <stdio.h>
#include <mpi.h>

void main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    int p;
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    printf("My rank is %d of %d processors\n", my_rank,p);
     
    MPI_Finalize();

}