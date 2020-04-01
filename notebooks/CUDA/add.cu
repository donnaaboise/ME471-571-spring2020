#include <stdio.h>

__device__ int addem( int a, int b ) 
{
    return a + b;
}

__global__ void add( int a, int b, int *c ) 
{
    *c = addem( a, b );
}

int main(void) 
{
    int a,b,c;
    int *dev_c;

    /* Allocate memory on the device */
    cudaMalloc( (void**)&dev_c, sizeof(int));

    a = 2;
    b = 7;
    add<<<1,1>>>(a, b, dev_c );

    /* Copy contens of dev_c back to c */
    cudaMemcpy( &c, dev_c, sizeof(int), cudaMemcpyDeviceToHost);

    printf( "%d + %d = %d\n", a,b,c );

    cudaFree( dev_c);

}


