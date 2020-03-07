
#include <mpi.h>
#include <math.h>

#include <stdio.h>
#include <stdlib.h>

void main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);


    /* ------------------------------------ User input ---------------------------------*/
    int N = atoi(argv[1]);
    int use_gatherv = atoi(argv[2]);

    double a = 0;
    double b = 1;

    /* ------------------------------- Numerical parameters ----------------------------*/

    double Nlocal = N/nprocs;
    double w = (b-a)/nprocs;

    double dx = (b-a)/N;
    double dx2 = dx*dx;

    /* Initialize data */
    double *const qmem = (double*) malloc((Nlocal+3)*sizeof(double));
    double *const q = &qmem[1];

    /* --------------------------------- Initialize data -------------------------------*/

    double pi = 4.*atan(1.0);
    for(int i = -1; i <= Nlocal+1; i++)
    {
        double x = a + rank*w + i*dx;
        q[i] = cos(2*pi*x);        
    }

    /* ----------------------------- Write data to a file  -----------------------------*/

    int node0 = 0;
    FILE *fout;
    double *qbig;
    int *recvcounts;
    int *displs;

    /* Write out meta-data */
    if (rank == 0)
    {
        fout = fopen("gather.out","w");
        fwrite(&N,1,sizeof(int),fout);
        fwrite(&a,1,sizeof(double),fout);
        fwrite(&b,1,sizeof(double),fout);
        fwrite(&use_gatherv,1,sizeof(int),fout);

        qbig = (double*) malloc((N+1)*sizeof(double));
        if (use_gatherv)
        {
            recvcounts = malloc(nprocs*sizeof(int));
            displs = malloc(nprocs*sizeof(int));

            for(int p = 0; p < nprocs; p++)
            {
                displs[p] = p*Nlocal;
                recvcounts[p] = Nlocal+1;
            }
        }
    }

    if (use_gatherv)
    {
        MPI_Gatherv(&q[0],Nlocal+1,MPI_DOUBLE,
                    qbig,recvcounts,displs,MPI_DOUBLE,
                    node0,MPI_COMM_WORLD);        
    }
    else
    {
        MPI_Gather(&q[0],Nlocal+1,MPI_DOUBLE,qbig,Nlocal+1,MPI_DOUBLE,node0,MPI_COMM_WORLD);        
    }

    if (rank == 0)
    {
        fwrite(qbig,N+1,sizeof(double),fout);
        free(qbig);
        fclose(fout);
    }
    free(qmem);
    
    MPI_Finalize();
}
