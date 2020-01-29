#include "demo4.h"
#include "demo_util.h"

#include <mpi.h>

void main(int argc, char** argv)
{
    /* Data arrays */
    double *x;
    double s, total_sum;
    double mean;

    /* MPI variables */
    int my_rank, nprocs;
    int tag = 0;
    int count;
    int source;
    MPI_Status status;

    /* Other local variables */
    int i, p;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    set_rank(my_rank);  /* Needed to log printing */

    /* Hardwire length of array */
    int n_global = pow2(10);
    int n_local = n_global/nprocs;

    if (my_rank == 0)
    {
        random_array(n_global,&x);  
        int p;
        for(p = 1; p < nprocs; p++)
        {
            int dest = p;
            MPI_Send((void*) &x[p*n_local],n_local,MPI_DOUBLE,dest,tag,MPI_COMM_WORLD);
        }
    }
    else
    {
        source = 0;
        empty_array(n_local,&x);
        MPI_Recv((void*) x,n_local,MPI_DOUBLE,source,tag,MPI_COMM_WORLD,&status);
        MPI_Get_count(&status,MPI_DOUBLE,&count);
        if (count < n_local)
        {
            print_debug("Something went wrong\n");
            print_debug("n_local = %d\n",n_local);
            print_debug("count = %d\n",count);
            exit(0);
        }
    }


    int root = 0;
    s = sum_array(n_local,x);   

/* Use #def 1 to get other branch of pragma statement */
#if 0
    print_global("Using MPI_Reduce + MPI_Bcast\n");
    MPI_Reduce(&s,&total_sum,1,MPI_DOUBLE,MPI_SUM,root,MPI_COMM_WORLD);    
    mean = total_sum/n_global;

    print_debug("The mean is  %.16f\n",my_rank,mean);        

    root = 0;
    MPI_Bcast(&mean,1,MPI_DOUBLE,root,MPI_COMM_WORLD);

    print_debug("The mean is %.16f\n",my_rank,mean);
#else
    print_global("Using MPI_Allreduce\n");
    MPI_Allreduce(&s,&total_sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);    
    mean = total_sum/n_global;

    print_debug("The mean is  %.16f\n",my_rank,mean);        
#endif

    if (my_rank == 0)
    {
        double mean_true;
        mean_true = sum_array(n_global,x)/n_global;
        print_global("True mean is %.16f\n",my_rank,mean_true);
    }

    /* Fill in details for broadcasting to all processors, computing 
    sums needed for STD, and then reducing results */

    delete_array(&x);

    MPI_Finalize();
}