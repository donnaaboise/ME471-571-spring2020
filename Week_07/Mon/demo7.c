#include <stdio.h>
#include <mpi.h>

#include "array_util.h"

void main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    int tag = 0;
    random_seed();


    MPI_Request request1;
    MPI_Request request2;
    if (my_rank == 0)
    {
        int p = 1;
        int dest = p;
        {
            double x = random_number();            
            printf("Sending a float %f to rank %d\n",x,p);
            MPI_Isend(&x,1,MPI_DOUBLE,dest,tag,MPI_COMM_WORLD,&request1);                
            printf("Done sending %f to rank %d\n",x,p);
        }
        {
            int a = 12345;
            printf("Sending an int %d to rank %d\n",a,p);
            MPI_Isend(&a,1,MPI_INTEGER,dest,tag,MPI_COMM_WORLD,&request2);
            printf("Done sending %d to rank %d\n",a,p);
        }
    }
    else if (my_rank == 1)
    {
        int sender = 0;
        {
            double x;
            printf("Rank %d is waiting to receive a float\n",my_rank);
            MPI_Irecv(&x,1,MPI_DOUBLE,sender,tag,MPI_COMM_WORLD,&request1); 
            printf("Rank %d received %f\n",my_rank,x);               
        }
        {
            int a;
            printf("Rank %d is waiting to receive an int\n",my_rank);
            MPI_Irecv(&a,1,MPI_INTEGER,sender,tag,MPI_COMM_WORLD,&request2); 
            printf("Rank %d received %d\n",my_rank,a);               
        }
    }

    /* Non-blocking tests */
    int flag;
    MPI_Test(&request1,&flag,MPI_STATUS_IGNORE);
    if (flag == 0)
    {
        printf("Request 1 is not complete\n");
    }
    else
    {
        printf("Request 1 is complete\n");
    }
    MPI_Test(&request2,&flag,MPI_STATUS_IGNORE);
    if (flag == 0)
    {
        printf("Request 2 is not complete\n");
    }
    else
    {
        printf("Request 2 is complete\n");
    }

    MPI_Finalize();
}