#include <mpi.h>

void main(int argc, char** argv)
{
    /* no MPI before this line */
    MPI_Init(&argc, &argv);

    /* Do something exciting! */

    MPI_Finalize();
    /* No more MPI commands */
}