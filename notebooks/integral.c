#include "array_util.h"
#include <math.h>
#include <mpi.h>
#include <stdio.h>

double integrand(double x)
{
    return pow((x-1),2)*exp(-pow(x,2));
}

int main(int argc, char** argv)
{
    /* MPI variables */
    MPI_Status status;

    MPI_Init(&argc, &argv);

    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* Total number of intervals */

    int N = atoi(argv[1]);    

    /* Interval over which to compute integral */
    double a = 0;
    double b = 1;
    double int_length = (b-a)/nprocs;

    /* Create interval for each processor */
    double intsub[2] = {a + my_rank*int_length, a + (my_rank+1)*int_length};

    int n_local = N/nprocs;
    double x0 = intsub[0];
    double x1 = intsub[1];
    double h = (b-a)/N;
    double integral = 0;
    for(int i = 0; i < n_local; i++)
    {
        double x = intsub[0] + h*i;
        integral += integrand(x)*h;
    }

    double total_sum;
    int root = 0;
    MPI_Reduce(&integral,&total_sum,1,MPI_DOUBLE,MPI_SUM,root,MPI_COMM_WORLD);

    if (my_rank == 0)
    {        
        double error = fabs(total_sum - 0.3041759198043616);
        //printf("%8d %24.16f %12.4e\n",N,total_sum, error);

        char fname[16];
        sprintf(fname,"integral_%02d.out",nprocs);
        FILE *fout = fopen(fname,"a");        
        fwrite(&N,1,sizeof(int),fout);
        fwrite(&error,1,sizeof(double),fout);
        fclose(fout);
    }

    /* Aggregate sum over each interval on root node */

    MPI_Finalize();
    return 0;
}