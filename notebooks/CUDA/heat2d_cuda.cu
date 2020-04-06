#include <stdio.h>
#include <stdlib.h>

#include <math.h>

#define IDX(i,j,cols,mbc)   i*cols + j + mbc  

__global__ void heat2d_update(int Nx, int Ny, int mbc, double dx, double dy,
                            double dt, double *qmem, double** qrows, 
                            double *qpmem, double** qprows);

__global__ void setup_arrays2d_cuda(int Nx, int Ny, int mbc, 
                                    double *qmem, double** qrows);


void allocate_2d(int N, int M, int mbc, double ***q)
{
    int rows = N + 1 + 2*mbc;
    int cols = M + 1 + 2*mbc; 

    double   *qmem = (double*) malloc(rows*cols*sizeof(double));
    double **qrows = (double**) malloc(rows*sizeof(double*));

    for(int i = 0; i < rows; i++)
    {
        qrows[i] = &qmem[cols*i + mbc];
    }    
    *q = &qrows[mbc];
}

void delete_2d(int mbc, double ***q)
{
    free(&(*q)[-mbc][-mbc]);
    free(&(*q)[-mbc]);
    *q = NULL;
}

/* Initial condition */
double init(double x, double y)
{
    return exp(-10*(x*x + y*y));
}

int main(int argc, char** argv)
{
    /* ------------------------------ Input parameters -------------------------------- */
    int Nx = atoi(argv[1]);
    int Ny = atoi(argv[2]);
    int nout = atoi(argv[3]);    

    /* Always print first and last time step, at a minimum */
    nout = (nout < 2) ? 2 : nout;  

    double L = 1;
    double Tfinal = 0.1;

    /* --------------------------- Numerical parameters ------------------------------- */
    double ax = -L;
    double bx = L;

    double ay = -L;
    double by = L;

    double dx = (bx-ax)/Nx;
    double dy = (by-ay)/Ny;
    double dx2 = dx*dx;
    double dy2 = dy*dy;

    /* Write out meta data  */
    FILE *fout = fopen("heat2d.out","w");        
    fwrite(&Nx,1,sizeof(int),fout);
    fwrite(&Ny,1,sizeof(int),fout);
    fwrite(&ax,1,sizeof(double),fout);
    fwrite(&bx,1,sizeof(double),fout);
    fwrite(&ay,1,sizeof(double),fout);
    fwrite(&by,1,sizeof(double),fout);

    /* ---------------------------- Initialize solution ------------------------------- */

    int mbc = 1;

    double **q;
    allocate_2d(Nx,Ny,mbc,&q);

    for(int i = -1; i <= Nx+1; i++)
    {
        double x = ax + i*dx;
        for(int j = -1; j <= Ny+1; j++)
        {
            double y = ay + j*dy;
            q[i][j] = init(x,y);            
        }
    }

    /* ----------------------------- Compute time step ---------------------------------*/
    /* Compute a stable time step
       1.  Estimate a stable time step 'dt_stable'.   This may not evenly divide Tfinal. 
       2.  Compute a minimum number M of time steps we need to take.
       3.  Divide Tfinal by M to get get a dt that is guaranteed smaller than dt_est and  
           satisfies M*dt = Tfinal.
    */
    double dsmin = (dx < dy) ? dx : dy;
    double ds2 = dsmin*dsmin;
    double dt_stable = 0.95*ds2/4;        /* Stable time step */
    int M = ceil(Tfinal/dt_stable) + 1;   /* Compute M to guarantee we hit Tfinal */
    double dt = Tfinal/M;                 /* dt <= dt_stable; M*dt = Tfinal */

    /* More meta data */
    fwrite(&M,1,sizeof(int),fout);

    /* 
       Set up an array of 'n' values that tell us when to save our solution 
       so that we save exactly nout time steps at roughly equally spaced times.  
    */
    int *noutsteps = (int*) malloc(nout*sizeof(int));    
    double dM = ((double) M-1)/(nout-1);
    dM = (dM < 1) ? 1 : dM;
    for(int m = 0; m <= nout-1; m++)
    {
        noutsteps[m] = (int) floor(m*dM);
    }

    /* Output initial condition */
    double t = 0;
    int k = 0;
    fwrite(&t,1,sizeof(double),fout);    
    for(int i = 0; i <= Nx; i++)
    {
        fwrite(&q[i][0],Ny+1,sizeof(double),fout);        
    }   
    k++;  /* Number of output files created */

    /* ---------------------------- Setup CUDA arrays ----------------------------------*/
    double *dev_q, *dev_qp;

    int rows  = (Nx + 1 + 2*mbc);
    int cols  = (Ny + 1 + 2*mbc);

    cudaMalloc( (void**) &dev_q, rows*cols*sizeof(double));
    double **dev_qrows;
    cudaMalloc( (void***) &dev_qrows, rows*sizeof(double*));
    setup_arrays2d_cuda<<<1,1>>>(Nx,Ny,mbc,dev_q, dev_qrows);

    cudaMalloc( (void**) &dev_qp, rows*cols*sizeof(double));
    double **dev_qprows;
    cudaMalloc( (void***) &dev_qprows, rows*sizeof(double*));
    setup_arrays2d_cuda<<<1,1>>>(Nx,Ny,mbc,dev_qp, dev_qprows);

    /* --------------------------- Start time stepping ---------------------------------*/
    /* Store q^{n+1} */
    double **qp;
    allocate_2d(Nx,Ny,mbc,&qp);


    /* Time loop;  compute q^{n+1} at each time step  */
    for(int n = 0; n <= M-1; n++)
    {
        t += dt;

        /* No-flux boundary conditions */
        for(int j = 0; j <= Ny; j++)
        {
            q[-1][j]   = q[1][j];
            q[Nx+1][j] = q[Nx-1][j];
        }
        for(int i = 0; i <= Nx; i++)
        {
            q[i][-1]   = q[i][1];
            q[i][Ny+1] = q[i][Ny-1];
        }

        /* ------------------------------ CUDA Kernel call -----------------------------*/

        int qsize = rows*cols;
        cudaMemcpy(dev_q, &q[-mbc][-mbc], qsize*sizeof(double), cudaMemcpyHostToDevice);


        int msize = Nx + 1;
        int nsize = Ny + 1;
        dim3 block(16,16);
        dim3 grid((msize+block.x - 1)/block.x, (nsize+block.y - 1)/block.y);

        heat2d_update<<<grid,block>>>(Nx, Ny, mbc, dx2, dy2, dt,
                                      dev_q, dev_qrows, 
                                      dev_qp, dev_qprows);

        cudaDeviceSynchronize();

        cudaMemcpy(&qp[-mbc][-mbc], dev_qp, qsize*sizeof(double), cudaMemcpyDeviceToHost);

        /* -------------------------------- Write output -------------------------------*/
        if (n == noutsteps[k])
        {            
            fwrite(&t,1,sizeof(double),fout);       
            for(int i = 0; i <= Nx; i++)
            {
                fwrite(&qp[i][0],Ny+1,sizeof(double),fout);        
            }   
            k++;
        }
        double **qtmp = q;
        q = qp;  
        qp = qtmp;
    }
    fclose(fout);

    cudaFree(dev_q);
    cudaFree(dev_qrows);

    cudaFree(dev_qp);
    cudaFree(dev_qprows);

    delete_2d(mbc,&q);
    delete_2d(mbc,&qp);
    free(noutsteps);

    return 0;
}