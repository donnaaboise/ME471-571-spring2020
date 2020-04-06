#include <math.h>
#include <stdlib.h>   /* For atoi */
#include <stdio.h>

#include "cg_serial.h"

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

void cg_1d(int N, int mbc, int kmax, double tol, int prt, matmult_t matmult, 
        double *F, double *u, int* it_cnt)
{

    /* ----------------------------------------------------------------
       Set up arrays and other vars needed for iterative method
    ---------------------------------------------------------------- */
    double *Auk  = allocate_1d(N,mbc);
    double *rk   = allocate_1d(N,mbc);
    double *rkp1 = allocate_1d(N,mbc);
    double *uk   = allocate_1d(N,mbc);
    double *ukp1 = allocate_1d(N,mbc);
    double *pk   = allocate_1d(N,mbc);
    double *wk   = allocate_1d(N,mbc);

    /* ----------------------------------------------------------------
       Start iterations
    ---------------------------------------------------------------- */
    for(int i = 0; i < N+1; i++)
    {
        uk[i] = 0;    /* or use memset */
        pk[i] = 0;
    }

    matmult(N,uk,Auk);

    for(int i = 1; i < N; i++)
    {
        /* Compute residual rk - F - A*uk */
        rk[i] = F[i] - Auk[i];
        pk[i] = rk[i];        
    }        

    *it_cnt = 0;
    double res;

    for(int k = 0; k < kmax; k++)
    {
        matmult(N,pk,wk);

        double a[2] = {0,0};
        for(int i = 1; i < N ; i++)
        {
            a[0] += rk[i]*rk[i];
            a[1] += pk[i]*wk[i];
        }
        double alpha = a[0]/a[1];

        double b[2] = {0,0};
        double norm_zk = 0, zk;
        for(int i = 1; i < N; i++)
        {
            double zk = alpha*pk[i];
            ukp1[i] = uk[i] + zk;
            rkp1[i] = rk[i] - alpha*wk[i];
            b[0] += rkp1[i]*rkp1[i];
            b[1] += rk[i]*rk[i];
            norm_zk = fabs(zk) > norm_zk ? fabs(zk) : norm_zk;
        }
        double beta = b[0]/b[1];

        if (prt != 0)
        {
            printf("%8d %16.8e\n",k,norm_zk);            
        }

        /* save results for output */
        *it_cnt = k+1;
        res = norm_zk;

        for(int i = 1; i < N; i++)
        {
            pk[i] = rkp1[i] + beta*pk[i];
            rk[i] = rkp1[i];
            uk[i] = ukp1[i];
        }
        if (norm_zk < tol)
        {
            break;
        }
    }

    for(int i = 1; i < N; i++)
    {
        u[i] = uk[i];
    }

    delete_1d(mbc,(double**) &Auk);
    delete_1d(mbc,(double**) &rk);
    delete_1d(mbc,(double**) &rkp1);
    delete_1d(mbc,(double**) &uk);
    delete_1d(mbc,(double**) &ukp1);
    delete_1d(mbc,(double**) &pk);
    delete_1d(mbc,(double**) &wk);
}

