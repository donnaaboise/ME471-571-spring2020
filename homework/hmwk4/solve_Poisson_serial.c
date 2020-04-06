#include <math.h>
#include <stdlib.h>   /* For atoi */
#include <stdio.h>

#include "cg_serial.h"   /* For f_t def */

#define PI 3.14159265358979323846264338327

typedef double (*f_t)(double x);


static
double* allocate_1d(int N, int mbc)
{
    int rows = N + 1 + 2*mbc;
    double   *qmem = malloc(rows*sizeof(double));
    return &qmem[mbc];
}

static
void delete_1d(int mbc, double **q)
{
    free(&(*q)[-mbc]);
    *q = NULL;
}


double utrue1(double x)
{
    double pi2 = 2*PI;
    double u = sin(pi2*x);
    return u;
}

double rhs1(double x)
{
    double pi2 = 2*PI;
    double fx = -(pi2)*(pi2)*sin(pi2*x);
    return fx;
}

double utrue2(double x)
{
    double pi = PI;
    double a = 40;

    double f = exp(-a*pow(x-0.5,2.0));
    double g = sin(2*pi*x);

    double u = f*g;
    return u;
}

double rhs2(double x)
{
    double pi = PI;
    double a = 40;

    double f = exp(-a*pow(x-0.5,2.0));
    double fp = -2*a*(x-0.5)*f;
    double fpp = -2*a*(x-0.5)*fp - 2*a*f;

    double g = sin(2*pi*x);
    double gp = 2*pi*cos(2*pi*x);
    double gpp = -pow(2*pi,2.0)*sin(2*pi*x);
    double u = f*g;
    double uxx = fpp*g + 2*fp*gp + f*gpp;
    return uxx;
}

static double s_dx;

void matmult(int N, double* x,double *Ax)
{
    double dx2 = s_dx*s_dx;
    x[0] = 0;
    x[N] = 0;
    for(int i = 1; i < N; i++)
    {
        Ax[i] = (x[i-1] - 2*x[i] + x[i+1])/dx2;            
    }    
}


void main(int argc, char** argv)
{
    /* ----------------------------------- user defined parameters -------------------- */
    int N = atoi(argv[1]);
    double tol = atof(argv[2]);
    double kmax = atoi(argv[3]);
    int prt = atoi(argv[4]);

    int mbc = 0;  /* Number of ghost cells */

    double a = 0;
    double b = 1;

    /* Assign true right hand side (u_xx) and solution */
    f_t utrue = utrue2;
    f_t rhs = rhs2;

    /* ----------------------------------- Numerical parameters ----------------------- */
    double dx = (b-a)/N;
    s_dx = dx;   /* Set static variable to be used in matmult */

    /* -------------------------------------------- RHS ------------------------------- */

    double *const F = allocate_1d(N,mbc);
    for(int i = 1; i < N; i++)
    {
        double x = a + i*dx;
        F[i] = rhs(x);
    }
    double dx2 = dx*dx;

    /* Include non-zero boundary conditions */
    F[1] -= utrue(a)/dx2;
    F[N-1] -= utrue(b)/dx2;

    /* ----------------------------------------------------------------
        Solve Poisson Problem using CG
    ---------------------------------------------------------------- */

    double *const u = allocate_1d(N,mbc);
    int it_cnt;
    cg_1d(N,mbc,kmax, tol, prt, matmult,&F[0],&u[0],&it_cnt);

    /* ----------------------------------------------------------------
       Calculate error and report results
    ---------------------------------------------------------------- */
    double err[3] = {0,0,0};
    for(int i = 1; i < N; i++)
    {
        double x = a + i*dx;
        double udiff = u[i] - utrue(x);
        err[0] += fabs(udiff)*dx;
        err[1] += fabs(udiff*udiff)*dx;
        err[2] = fabs(udiff) > err[2] ? fabs(udiff) : err[2];
    }
    err[1] = sqrt(err[1]);    /* 2-norm */

    printf("%10d %10d %12.4e %12.4e %12.4e\n",N,it_cnt, 
                 err[0],err[1],err[2]);

    delete_1d(mbc,(double**) &F);
    delete_1d(mbc,(double**) &u);
}

