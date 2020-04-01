#include <stdio.h>

__global__ void add( int *c) 
{
    int id = blockIdx.x;
    c[id] = id;
}

#define N 10

int main(void) 
{
    int c[N];
    int *dev_c;
    int i;

    /* Allocate memory on the device */
    cudaMalloc( (void**)&dev_c, N*sizeof(int));

    add<<<N,1>>>(dev_c );

    /* Copy contents of dev_c back to c */
    cudaMemcpy( &c, dev_c, N*sizeof(int), cudaMemcpyDeviceToHost);

    for(i = 0; i < N; i++)
    {
        printf( "c[%d] = %d\n",i,c[i]);
    }

    cudaFree(dev_c);

}


