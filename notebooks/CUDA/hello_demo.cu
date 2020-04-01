#include <stdio.h>

__global__ void kernel( void ) 
{
    /* Do something fun! */
}

int main(void) 
{
    kernel<<<1,1>>>();

    printf( "Hello, World!\n" );

    return 0;
}
