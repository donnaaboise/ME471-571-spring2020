#include <mpi.h>
#include <math.h>

#include <stdio.h>
#include <stdlib.h>

typedef double (*f_t)(double x);

/* Function to differentiate */
double f1(double x)
{
    double pi = 4.*atan(1.0);
    return sin(2*pi*x);
}

/* Function to differentiate */
double f1_deriv1(double x)
{
    double pi = 4.*atan(1.0);
    return 2*pi*cos(2*pi*x);
}

/* Function to differentiate */
double f1_deriv2(double x)
{
    double pi = 4.*atan(1.0);
    return -pow(2*pi,2.)*sin(2*pi*x);
}


/* Integrate this over [0,1] */
double f2(double x)
{
    return 1/(1 + 25*x*x);
}

void main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);


    /* ------------------------------------ User input ---------------------------------*/
    int N = atoi(argv[1]);
    int deriv_choice = atoi(argv[2]);

    double a = 0;
    double b = 1;

    f_t f = f1;  /* Function pointer */
    f_t fderiv;
    switch(deriv_choice)
    {
        case 1: 
            fderiv = f1_deriv1;
            break;
        case 2:
            fderiv = f1_deriv2;
            break;            
        default:
            printf("Derivatives > 2 not yet supported.\n");
            exit(0);
    }

    /* ------------------------------- Numerical parameters ----------------------------*/

    double m = N/nprocs;
    double w = (b-a)/nprocs;

    double dx = (b-a)/N;
    double dx2 = dx*dx;

    /* Initialize data */
    double *qmem = (double*) malloc((m+3)*sizeof(double));
    double *smem = (double*) malloc((m+3)*sizeof(double));

    double *const q = &qmem[1];
    double *const soln = &smem[1];

    for(int i = -1; i <= m+1; i++)
    {
        double x = a + rank*w + i*dx;
        q[i] = f(x);
        soln[i] = fderiv(x);
    }

    /* ----------------------------- Write data to a file  -----------------------------*/

    int node0 = 0;
    FILE *fout;
    double *qbig;
    int *recvcounts;
    int *displs;

    /* Write out some meta-data and the true derivative */
    if (rank == 0)
    {
        fout = fopen("deriv.out","w");
        fwrite(&N,1,sizeof(int),fout);
        fwrite(&a,1,sizeof(double),fout);
        fwrite(&b,1,sizeof(double),fout);
        fwrite(&deriv_choice,1,sizeof(int),fout);

        qbig = (double*) malloc((N+1)*sizeof(double));
        recvcounts = malloc(nprocs*sizeof(int));
        displs = malloc(nprocs*sizeof(int));

        for(int p = 0; p < nprocs; p++)
        {
            displs[p] = p*m;
            recvcounts[p] = m+1;
        }
    }

#if 0
    MPI_Gatherv(&soln[0],m+1,MPI_DOUBLE,
                qbig,recvcounts,displs,MPI_DOUBLE,
                node0,MPI_COMM_WORLD);
#else    
    MPI_Gather(&soln[1],m,MPI_DOUBLE,qbig,m,MPI_DOUBLE,node0,MPI_COMM_WORLD);
#endif

    if (rank == 0)
    {
        fwrite(qbig,N+1,sizeof(double),fout);
    }

    /* --------------------------- Compute the derivative  -----------------------------*/

    double *dmem = (double*) malloc((N+3)*sizeof(double));  
    double *deriv = &dmem[1];

    /* Compute approximate derivatives */
    for(int i = 0; i <= N; i++)
    {
        /* first derivative */
        switch(deriv_choice)
        {
            case 1: 
                /* First Derivative */
                deriv[i] = (q[i+1] - q[i-1])/(2*dx);
                break;
            case 2:
                /* Second derivative */
                deriv[i] = (q[i+1] - 2*q[i] +  q[i-1])/dx2;
                break;
            default:
                printf("Only derivatives 1 and 2 supported at this time.\n");
                exit(0);               
        }
    }

#if 0
    MPI_Gatherv(&deriv[0],m+1,MPI_DOUBLE,
                qbig,recvcounts,displs,MPI_DOUBLE,
                node0,MPI_COMM_WORLD);
#else                
    MPI_Gather(&deriv[0],m,MPI_DOUBLE,qbig, m,MPI_DOUBLE,node0,MPI_COMM_WORLD);
#endif    

    if (rank ==0)
    {
        fwrite(qbig,N+1,sizeof(double),fout);        
        free(qbig);
        fclose(fout);
    }

    free(qmem);
    free(dmem);
    free(smem);

    MPI_Finalize();
}