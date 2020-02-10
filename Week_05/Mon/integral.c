#include "array_util.h"
#include <math.h>
#include <mpi.h>
#include <stdio.h>

double integrand(double x)
{
    return pow((x-1),2)*exp(-pow(x,2));
}

void main(int argc, char** argv)
{
    /* MPI variables */
    MPI_Status status;

    MPI_Init(&argc, &argv);

    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    set_rank(my_rank);  /* Needed to log printing */

    /* Total number of intervals */
    int N = pow2(10);

    /* Interval over which to compute integral */
    double a = 0;
    double b = 1;
    double int_length = (b-a)/nprocs;

    int tag = 0;
    double intsub[2];
    if (my_rank == 0)
    {
        /* Compute integral */
        for(int p = 1; p < nprocs; p++)
        {
            int dest = p;
            intsub[0] = a + p*int_length;
            intsub[1] = a + (p+1)*int_length;
            MPI_Send(intsub,2,MPI_DOUBLE,dest,tag,MPI_COMM_WORLD);
        }
    }
    else
    {
        int source = 0;
        MPI_Recv(intsub,2,MPI_DOUBLE,source,tag,MPI_COMM_WORLD,&status);
        /* Receive interval information from each processor */
    }

    /* Create array using interval */

    /* solution : 0.3041759198043616 */
    int n_local = N/nprocs;
    double integral = 0;
    double h = (b-a)/N;
    for(int i = 0; i < n_local; i++)
    {
        double x = intsub[0] + h*i;
        integral += integrand(x);
    }

    int dest = 0;
    MPI_Send(&integral,1,MPI_DOUBLE,dest,tag,MPI_COMM_WORLD);

    if (my_rank == 0)
    {
        /* Compute integral */
        double total_sum = integral;
        for(int p = 0; p < nprocs; p++)
        {
            int source = p;
            double integral;
            MPI_Recv(&integral,1,MPI_DOUBLE,source,tag,MPI_COMM_WORLD,&status);
            total_sum += integral;
        }
    }
    if (my_rank == 0)
    {        
        double error = total_sum - 0.3041759198043616;
        printf("%8d %24.16f %12.4e\n",N,total_sum, error);
    }

    /* Aggregate sum over each interval on root node */

    MPI_Finalize();
}