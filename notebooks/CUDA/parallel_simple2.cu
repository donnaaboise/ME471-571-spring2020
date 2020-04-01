#include <stdio.h>

__global__ void add( int *c) 
{
    int id = threadIdx.x;
    c[id] = id;
}

#define N 16

int main(void) 
{
    int c[N];
    int *dev_c;
    int i;

    /* Allocate memory on the device */
    cudaMalloc( (void**)&dev_c, N*sizeof(int));

    add<<<1,N>>>(dev_c);

    /* Copy contents of dev_c back to c */
    cudaMemcpy( &c, dev_c, N*sizeof(int), cudaMemcpyDeviceToHost);

    for(i = 0; i < N; i++)
    {
        printf( "c[%d] = %d\n",i,c[i]);
    }

    cudaFree(dev_c);

}


