#include <math.h>
#include <stdlib.h>   /* For atoi */
#include <stdio.h>
#include <mpi.h>
#include <string.h>

#include "solver.h"


void splitting(int N_global, int mbc, int kmax, double tol, int prt, 
               matmult2d_t matmult, method_t method, double **F, 
               double **u, int* it_cnt, double *res, MPI_Comm comm_cart)
{

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int maxdims = 2;
    int dims[2],periods[2], mycoords[2];
    MPI_Cart_get(comm_cart,maxdims,dims,periods,mycoords);

    int Nx = N_global/dims[0];
    int Ny = N_global/dims[1];

    double h = 1.0/N_global;
    double h2 = h*h;

    /* ----------------------------------------------------------------
       Set up arrays and other vars needed for iterative method
    ---------------------------------------------------------------- */
    double **Auk   = allocate_2d(Nx,Ny,mbc);
    double **uk    = allocate_2d(Nx,Ny,mbc);
    double **ukp1  = allocate_2d(Nx,Ny,mbc);
    double **rk    = allocate_2d(Nx,Ny,mbc);

    /* ----------------------------------------------------------------
       Start iterations
    ---------------------------------------------------------------- */

    int nsize = (Nx + 1 + 2*mbc)*(Ny + 1 + 2*mbc)*sizeof(double);
    memset(&uk[-mbc][-mbc],0,nsize);
    memset(&rk[-mbc][-mbc],0,nsize);
    memset(&Auk[-mbc][-mbc],0,nsize);
    memset(&ukp1[-mbc][-mbc],0,nsize);

    *it_cnt = 0;

    for(int k = 0; k < kmax; k++)
    {
        matmult(N_global,uk,Auk, comm_cart);

        for(int i = 0; i < Nx+1; i++)
        {
            for(int j = 0; j < Ny+1; j++)
            {
                if (method == JA)
                {
                    rk[i][j] = F[i][j] - Auk[i][j];                    
                }
                else if (method == GS)
                {
                    rk[i][j] = F[i][j] - 
                        (ukp1[i-1][j] + uk[i+1][j] + ukp1[i][j-1] + uk[i][j+1] - 4*uk[i][j])/h2;
                }
                ukp1[i][j] = uk[i][j] - h2*rk[i][j]/4.0;
            }
        }

        double rmax = 0;
        for(int i = 0; i < Nx+1; i++)
        {
            for(int j = 0; j < Ny+1; j++)
            {
                double r = fabs(rk[i][j]);
                rmax = r > rmax ? r : rmax;
            }
        }

        double norm_rk;
        MPI_Allreduce(&rmax,&norm_rk,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);        

        if (prt && rank == 0)
        {
            printf("%5d %24.16e\n",k,norm_rk);
        }

        /* save results for output */
        *it_cnt = k+1;
        *res = norm_rk;

        if (norm_rk < tol)
        {
            break;
        }
        double **tmp;
        tmp = uk;
        uk = ukp1;
        ukp1 = tmp;
    }

    for(int i = 0; i < Nx+1; i++)
    {
        for(int j = 0; j < Ny+1; j++)
        {
            u[i][j] = uk[i][j];            /* Or just use memcpy */
        }
    }

    delete_2d(mbc,(double***) &Auk);
    delete_2d(mbc,(double***) &uk);
    delete_2d(mbc,(double***) &ukp1);
}

