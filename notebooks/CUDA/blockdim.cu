#include <stdio.h>
#include <sys/time.h>


#define CHECK(call)                                                          \
{                                                                            \
    const cudaError_t error = call;                                          \
    if (error != cudaSuccess)                                                \
    {                                                                        \
        printf("Error : %s : %d, ", __FILE__, __LINE__);                     \
        printf("code : %d, reason : %s\n",error, cudaGetErrorString(error)); \
        exit(1);                                                             \
    }                                                                        \
}       

void allocate_2d(int N, int M, double ***q)
{
    int rows = N;
    int cols = M;

    double   *qmem = (double*) malloc(rows*cols*sizeof(double));
    double **qrows = (double**) malloc(rows*sizeof(double*));

    for(int i = 0; i < rows; i++)
    {
        qrows[i] = &qmem[cols*i];
    }    
    *q = &qrows[0];
}

void delete_2d(double ***q)
{
    free(&(*q)[0][0]);
    free(&(*q)[0]);
    *q = NULL;
}


__global__ void setup_arrays2d_cuda(int Nx, int Ny,
                                    double *qmem, double** qrows, double*** q)
{
    int rows = Nx;
    int cols = Ny; 

    for(int i = 0; i < rows; i++)
    {
        qrows[i] = &qmem[cols*i];
    }    
    *q = &qrows[0];
}


__global__ void copymat_x(int Nx, int Ny, double*** dev_A, double ***dev_B) 
{    
    double **A = *dev_A;
    double **B = *dev_B;

    /* This thread runs over i values (rows) */
    int j = threadIdx.x + blockIdx.x*blockDim.x;
    if (j < Ny) 
        for(int i = 0; i < Nx; i++) 
            B[i][j]   = A[i][j];
}

__global__ void copymat_y(int Nx, int Ny, double*** dev_A, double*** dev_B) 
{    
    double **A = *dev_A;
    double **B = *dev_B;

    /* This thread runs over j values (columns) */
    int i = threadIdx.x + blockIdx.x*blockDim.x;
    if (i < Nx)
        for(int j = 0; j < Ny; j++) 
            B[i][j]   = A[i][j];        
}

double cpuSecond()
{
    struct timeval tp;
    gettimeofday(&tp,NULL);
    return (double) tp.tv_sec + (double)tp.tv_usec*1e-6;
}

int main(int argc, char* argv[]) 
{

    int P = atoi(argv[1]);
    int block_dim = atoi(argv[2]);

    if (P > 13)
        printf("P is probably too large and will likely seg. fault\n");

    int Nx = 1 << P;  
    int Ny = 1 << P;  

    size_t nbytes = Nx*Ny*sizeof(double);
    printf("N = %d\n",Nx);
    printf("Memory = %dMb\n",nbytes/(1024*1024));

    double *dev_Amem, **dev_Arows, ***dev_A;
    cudaMalloc( (void**) &dev_Amem, Nx*Ny*sizeof(double));
    cudaMalloc( (void***) &dev_Arows,  Ny*sizeof(double*));
    cudaMalloc( (void****) &dev_A,        sizeof(double**));
    setup_arrays2d_cuda<<<1,1>>>(Nx,Ny,dev_Amem, dev_Arows,dev_A);    

    double *dev_Bmem, **dev_Brows, ***dev_B;
    cudaMalloc( (void**) &dev_Bmem, Nx*Ny*sizeof(double));
    cudaMalloc( (void***) &dev_Brows,  Ny*sizeof(double*));
    cudaMalloc( (void****) &dev_B,        sizeof(double**));
    setup_arrays2d_cuda<<<1,1>>>(Nx,Ny,dev_Bmem, dev_Brows,dev_B);    


    double **A;
    allocate_2d(Nx, Ny, &A);
    memset(&A[0][0],0,nbytes);

    cudaMemcpy(dev_Amem, &A[0][0], nbytes, cudaMemcpyHostToDevice);

    dim3 block(block_dim);  
    dim3 grid((Ny+block.x-1)/block.x);
    printf("Number of blocks : %d\n",grid.x);

    double etime[2];

    /* Run over rows */
    double start = cpuSecond();
    copymat_x<<<grid,block>>>(Nx,Ny,dev_A, dev_B);
    CHECK(cudaDeviceSynchronize());
    etime[0] = cpuSecond() - start;
    printf("GPU Kernel (rows) %12.6f (s)\n",etime[0]);

    start = cpuSecond();
    copymat_y<<<grid,block>>>(Nx,Ny,dev_A, dev_B);
    CHECK(cudaDeviceSynchronize());
    etime[1] = cpuSecond() - start;
    printf("GPU Kernel (cols) %12.6f (s)\n",etime[1]);

    printf("Ratio (cols/rows) %12.2f\n",etime[1]/etime[0]);

    cudaFree(dev_Amem);
    cudaFree(dev_Arows);
    cudaFree(dev_A);

    cudaFree(dev_Bmem);
    cudaFree(dev_Brows);
    cudaFree(dev_B);

    delete_2d(&A);

    cudaDeviceReset();
}


