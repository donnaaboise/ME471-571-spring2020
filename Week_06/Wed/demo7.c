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
    if (my_rank == 0)
    {
        /* Send a message to each proc! */
        for(int p = 1; p < nprocs; p++)
        {
            int dest = p;
            {
                double x = random_number();            
                printf("Sending %f to rank %d\n",x,p);
                MPI_Isend(&x,1,MPI_DOUBLE,dest,tag,MPI_COMM_WORLD,&request1);                
                printf("Done sending %f to rank %d\n",x,p);
            }
            {
                int a = 12345;
                printf("Sending %d to rank %d\n",a,p);
                MPI_Isend(&a,1,MPI_INTEGER,dest,tag,MPI_COMM_WORLD,&request2);
                printf("Done sending %d to rank %d\n",a,p);
            }

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
        }
    }

    


    MPI_Finalize();
}