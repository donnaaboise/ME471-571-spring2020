#include <stdio.h>
#include <mpi.h>
#include <time.h>

void sleep(double t_total)
{
    double t0, t1;
    t0 = clock();
    t1 = t0;
    while ((t1-t0)/CLOCKS_PER_SEC < t_total)
    {
        t1 = clock();
    }
}

void main(int argc, char** argv)
{
    int my_rank;
    int p;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    MPI_Comm_size(MPI_COMM_WORLD, &p);

    printf("My rank is %d\n", my_rank);
     
    sleep(2.0);
    
    printf("Rank %d is done!\n",my_rank);    
    

    MPI_Finalize();

}