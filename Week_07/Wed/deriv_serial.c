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

    double dx = (b-a)/N;
    double dx2 = dx*dx;

    /* Initialize data */
    double *xmem = (double*) malloc((N+3)*sizeof(double));
    double *qmem = (double*) malloc((N+3)*sizeof(double));
    double *smem = (double*) malloc((N+3)*sizeof(double));

    double *x = &xmem[1];
    double *q = &qmem[1];
    double *soln = &smem[1];

    for(int i = -1; i < N+2; i++)
    {
        x[i] = a + i*dx;
        q[i] = f(x[i]);
        soln[i] = fderiv(x[i]);
    }

    /* ----------------------------- Write data to a file  -----------------------------*/

    FILE *fout = fopen("deriv.out","w");
    fwrite(&N,1,sizeof(int),fout);
    fwrite(&a,1,sizeof(double),fout);
    fwrite(&b,1,sizeof(double),fout);
    fwrite(&deriv_choice,1,sizeof(int),fout);

    fwrite(&soln[0],N+1,sizeof(double),fout);

    /* --------------------------- Compute the derivative  -----------------------------*/

    double *dmem = (double*) malloc((N+3)*sizeof(double));
    double *deriv = &dmem[1];

    for(int i = 0; i < N; i++)
    {
        /* first derivative */
        switch(deriv_choice)
        {
            case 1: 
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

    fwrite(&deriv[0],N+1,sizeof(double),fout);

    fclose(fout);

    free(xmem);
    free(qmem);
    free(dmem);
    free(smem);

    MPI_Finalize();
}