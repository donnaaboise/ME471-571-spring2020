#include <math.h>
#include <stdlib.h>   /* For atoi */
#include <stdio.h>
#include <string.h>
#include <mpi.h>

#include "solver.h"


static
double dot_product(int Nx, int Ny, double **u, double **v)
{
    double d = 0;
    for(int i = 0; i < Nx+1 ; i++)
        for(int j = 0; j < Ny+1; j++)
            d += u[i][j]*v[i][j];

    return d;
}

static
double vec_max(int Nx, int Ny, double **u)
{
    double dmax = 0;
    for(int i = 0; i < Nx+1 ; i++)
        for(int j = 0; j < Ny+1; j++)
            dmax =  fabs(u[i][j]) > dmax ? fabs(u[i][j]) : dmax;

    return dmax;
}


void cg_2d(int N_global,  int mbc, int kmax, double tol, int prt, matmult2d_t matmult, 
        double **F, double **u, int* it_cnt, double *res, MPI_Comm comm_cart)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int maxdims = 2;
    int dims[2],periods[2], mycoords[2];
    MPI_Cart_get(comm_cart,maxdims,dims,periods,mycoords);

    int Nx = N_global/dims[0];
    int Ny = N_global/dims[1];    

    /* ----------------------------------------------------------------
       Set up arrays and other vars needed for iterative method
    ---------------------------------------------------------------- */
    double **Auk  = allocate_2d(Nx,Ny,mbc);
    double **rk   = allocate_2d(Nx,Ny,mbc);
    double **rkp1 = allocate_2d(Nx,Ny,mbc);
    double **uk   = allocate_2d(Nx,Ny,mbc);
    double **ukp1 = allocate_2d(Nx,Ny,mbc);
    double **pk   = allocate_2d(Nx,Ny,mbc);
    double **zk   = allocate_2d(Nx,Ny,mbc);
    double **Apk  = allocate_2d(Nx,Ny,mbc);

    /* ----------------------------------------------------------------
       Start iterations
    ---------------------------------------------------------------- */

    int nsize = (Nx + 1 + 2*mbc)*(Ny + 1 + 2*mbc)*sizeof(double);
    memset( &Auk[-mbc][-mbc],0,nsize);
    memset(  &rk[-mbc][-mbc],0,nsize);
    memset(&rkp1[-mbc][-mbc],0,nsize);
    memset(  &uk[-mbc][-mbc],0,nsize);
    memset(&ukp1[-mbc][-mbc],0,nsize);    
    memset(  &pk[-mbc][-mbc],0,nsize);
    memset(  &zk[-mbc][-mbc],0,nsize);
    memset( &Apk[-mbc][-mbc],0,nsize);

    int i1, i2, j1, j2;
    get_endpoints_indices(N_global, &i1, &i2, &j1, &j2,comm_cart);

    matmult(N_global,uk,Auk,comm_cart);

    for(int i = 0; i < Nx+1; i++)
        for(int j = 0; j < Ny+1; j++)
        {
            /* Compute residual rk - F - A*uk */
            rk[i][j] = F[i][j] - Auk[i][j];
            pk[i][j] = rk[i][j];        
        }

    *it_cnt = 0;
    for(int k = 0; k < kmax; k++)
    {
        matmult(N_global,pk,Apk,comm_cart);

        double a_local[2];
        a_local[0] = dot_product(Nx,Ny, rk,rk);
        a_local[1] = dot_product(Nx,Ny, pk,Apk);
        double a[2];
        MPI_Allreduce(&a_local,&a,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);        

        double alpha = a[0]/a[1];

        for(int i = 0; i < Nx+1; i++)
            for(int j = 0; j < Ny+1; j++)
            {
                zk[i][j]   = alpha*pk[i][j];
                ukp1[i][j] = uk[i][j] + zk[i][j];
                rkp1[i][j] = rk[i][j] - alpha*Apk[i][j];
            }

        double dmax = vec_max(Nx,Ny,zk);
        double norm_zk;
        MPI_Allreduce(&dmax,&norm_zk,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);        

        double b_local[2];
        b_local[0] = dot_product(Nx,Ny,rkp1,rkp1);
        b_local[1] = a_local[0];
        double b[2];
        MPI_Allreduce(&b_local[0],&b[0],2,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);        
        double beta = b[0]/b[1];

        if (prt != 0 && rank == 0)
            printf("%8d %16.8e\n",k,norm_zk);            

        /* save results for output */
        *it_cnt = k+1;
        *res = norm_zk;

        /* Update the A-conjugate direction vectors pk */
        for(int i = 0; i < Nx+1; i++) 
            for(int j = 0; j < Ny+1; j++) 
                pk[i][j] = rkp1[i][j] + beta*pk[i][j];

        if (norm_zk < tol)
            break;

        double **tmp;
        tmp = rk;
        rk = rkp1;
        rkp1 = tmp;

        tmp = uk;
        uk = ukp1;
        ukp1 = tmp;
    }

    for(int i = 0; i < Nx+1; i++)
        for(int j = 0; j < Ny+1; j++)
            u[i][j] = uk[i][j];            /* Or just use memcpy */

    delete_2d(mbc,(double***) &Auk);
    delete_2d(mbc,(double***) &rk);
    delete_2d(mbc,(double***) &rkp1);
    delete_2d(mbc,(double***) &uk);
    delete_2d(mbc,(double***) &ukp1);
    delete_2d(mbc,(double***) &pk);
    delete_2d(mbc,(double***) &Apk);
}


