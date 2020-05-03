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
    {
        for(int j = 0; j < Ny+1; j++)
        {
            d += u[i][j]*v[i][j];
        }
    }    
    return d;
}

static
double vec_max(int Nx, int Ny, double **u)
{
    double dmax = 0;
    for(int i = 0; i < Nx+1 ; i++)
    {
        for(int j = 0; j < Ny+1; j++)
        {            
            dmax =  fabs(u[i][j]) > dmax ? fabs(u[i][j]) : dmax;
        }
    }    
    return dmax;
}


void bicg_2d(int N_global,  int mbc, int kmax, double tol, int prt, matmult2d_t matmult, 
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
    double **Apk  = allocate_2d(Nx,Ny,mbc);

    double **rk0   = allocate_2d(Nx,Ny,mbc);
    double **rhat  = allocate_2d(Nx,Ny,mbc);
    double **Arhat  = allocate_2d(Nx,Ny,mbc);

    /* ----------------------------------------------------------------
       Start iterations
    ---------------------------------------------------------------- */

    int nsize = (Nx + 1 + 2*mbc)*(Ny + 1 + 2*mbc)*sizeof(double);
    memset(&Auk[-mbc][-mbc],0,nsize);
    memset(&rk[-mbc][-mbc],0,nsize);
    memset(&rkp1[-mbc][-mbc],0,nsize);
    memset(&uk[-mbc][-mbc],0,nsize);
    memset(&ukp1[-mbc][-mbc],0,nsize);    
    memset(&pk[-mbc][-mbc],0,nsize);
    memset(&Apk[-mbc][-mbc],0,nsize);

    memset(&rk0[-mbc][-mbc],0,nsize);
    memset(&rhat[-mbc][-mbc],0,nsize);    
    memset(&Arhat[-mbc][-mbc],0,nsize);


    matmult(N_global,uk,Auk,comm_cart);

    for(int i = 0; i < Nx+1; i++)
    {
        for(int j = 0; j < Ny+1; j++)
        {
            /* Compute residual rk - F - A*uk */
            rk[i][j] = F[i][j] - Auk[i][j];
            pk[i][j] = rk[i][j];        
            rk0[i][j] = rk[i][j];
        }
    }        

    *it_cnt = 0;
    for(int k = 0; k < kmax; k++)
    {
        matmult(N_global,pk,Apk,comm_cart);

        /* Compute alpha */
        double a_local[2];
        a_local[0] = dot_product(Nx,Ny, rk,rk0);
        a_local[1] = dot_product(Nx,Ny, Apk,rk0);
        double a[2];
        MPI_Allreduce(a_local,a,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);   

        if (a[1] == 0)
            printf("rank %d : a[1] == 0\n",rank);

        double alpha = a[0]/a[1];

        /* Compute rhat */
        for(int i = 0; i < Nx+1; i++)
        {
            for(int j = 0; j < Ny+1; j++)
            {
                rhat[i][j] = rk[i][j] - alpha*Apk[i][j];
            }
        }
        matmult(N_global,rhat,Arhat,comm_cart);

        /* Compute omega */
        double w_local[2], w[2];
        w_local[0] = dot_product(Nx,Ny,rhat,Arhat);
        w_local[1] = dot_product(Nx,Ny,Arhat,Arhat);
        MPI_Allreduce(w_local,w,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);        
        if (w[1] == 0)
            printf("rank %d : w[1] == 0\n",rank);


        double omega = w[0]/w[1];

        /* Update solution and residual */
        double rmax = 0;
        for(int i = 0; i < Nx+1; i++)
        {
            for(int j = 0; j < Ny+1; j++)
            {
                uk[i][j] = uk[i][j] + alpha*pk[i][j] + omega*rhat[i][j];
                rkp1[i][j] = rhat[i][j] - omega*Arhat[i][j];
                rmax = (fabs(rkp1[i][j]) > rmax) ? fabs(rkp1[i][j]) : rmax;
            }
        }


        /* Compute beta */
        double b_local[2], b[2];
        b_local[0] = dot_product(Nx,Ny,rkp1,rk0);
        MPI_Allreduce(b_local,b,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);        
        b[1] = a[0];
        if (b[1] == 0)
            printf("rank %d : b[1] == 0\n",rank);

        double beta = (alpha/omega)*(b[0]/b[1]);


        /* Update directions pk */
        for(int i = 0; i < Nx+1; i++)
        {
            for(int j = 0; j < Ny+1; j++)
            {
                pk[i][j] = rkp1[i][j] + beta*(pk[i][j] - omega*Apk[i][j]);
            }
        }

        double norm_rkp1;
        MPI_Allreduce(&rmax,&norm_rkp1,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);        

        if (prt != 0 && rank == 0)
        {
            printf("%8d %16.8e\n",k,norm_rkp1);            
        }

        /* save results for output */
        *it_cnt = k+1;
        *res = norm_rkp1;

        if (norm_rkp1 < tol)
        {
            break;
        }

        double **tmp;
        tmp = rk;
        rk = rkp1;
        rkp1 = tmp;

#if 0
        tmp = uk;
        uk = ukp1;
        ukp1 = tmp;
#endif        
    }

    for(int i = 0; i < Nx+1; i++)
    {
        for(int j = 0; j < Ny+1; j++)
        {
            u[i][j] = uk[i][j];            /* Or just use memcpy */
        }
    }

    delete_2d(mbc,(double***) &Auk);
    delete_2d(mbc,(double***) &rk);
    delete_2d(mbc,(double***) &rkp1);
    delete_2d(mbc,(double***) &uk);
    delete_2d(mbc,(double***) &ukp1);
    delete_2d(mbc,(double***) &pk);
    delete_2d(mbc,(double***) &Apk);

    delete_2d(mbc,(double***) &rk0);
    delete_2d(mbc,(double***) &rhat);
    delete_2d(mbc,(double***) &Arhat);
}


