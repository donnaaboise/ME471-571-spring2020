#include <stdio.h>
#include <sys/time.h>

#define P (1<<16)

__global__ void copymat_x(int m, int n, int* A, int *B) 
{    
    int idx, ix;
    int iy = threadIdx.y + blockIdx.y*blockDim.y;
    if (iy < n) 
        for(ix = 0; ix < P; ix++) {
            idx  = iy*m + ix; 
            B[idx]   = A[idx];
        }
}

__global__ void copymat_y(int m, int n, int* A, int *B) 
{    
    int ix = threadIdx.x + blockIdx.x*blockDim.x;
    int idx, iy;
    if (ix < m)
        for(iy = 0; iy < P; iy++) {
            idx  = iy*m + ix; 
            B[idx]   = A[idx];        
        }
}

double cpuSecond()
{
    struct timeval tp;
    gettimeofday(&tp,NULL);
    return (double) tp.tv_sec + (double)tp.tv_usec*1e-6;
}

int main(int argc, char** argv) 
{
    size_t m = 1 << 16;  
    size_t n = 1 << 16;  
    size_t nbytes = m*n*sizeof(int);

    printf("P = %d\n",P);

    int* A = (int*) malloc(nbytes);
    int *B = (int*) malloc(nbytes);

    memset(A,0,nbytes);

    int *dev_A, *dev_B;
    cudaMalloc((void**) &dev_A, nbytes);
    cudaMalloc((void**) &dev_B, nbytes);
    cudaMemcpy(dev_A, A, nbytes, cudaMemcpyHostToDevice);

#if 0
    /* One thread per row */
    dim3 block(1,32);  
    dim3 grid(1,(n+block.y-1)/block.y);
    double start = cpuSecond();
    copymat_x<<<grid,block>>>(m,n,dev_A, dev_B);
#else
    /* One thread per column */
    dim3 block(32,1);  
    dim3 grid((m+block.x-1)/block.x,1);
    double start = cpuSecond();
    copymat_y<<<grid,block>>>(m,n,dev_A, dev_B);
#endif
    cudaDeviceSynchronize();
    double etime = cpuSecond() - start;
    printf("GPU Kernel %10.3g (s)\n",etime);

    cudaFree(dev_A);
    cudaFree(dev_B);
    free(A);
    free(B);

    cudaDeviceReset();
}


