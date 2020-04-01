#include <stdio.h>

int main( void ) 
{
    cudaDeviceProp  prop;

    int count;
    char str[4];
    cudaGetDeviceCount( &count);

    if (count == 0)
    {
        printf("No CUDA capable devices found.\n");
    }

    for (int i = 0; i < count; i++) 
    {
        cudaGetDeviceProperties( &prop, i);

        printf("   --- General Information for device %d ---\n", i );
        printf("Name:  %s\n", prop.name );
        printf("\n");
        sprintf(str,"%d.%d",prop.major, prop.minor);
        printf("Compute capability    :  %14s\n",str);
        printf("Clock rate            :  %14.2f (GHz)\n", prop.clockRate/1000000.0);
        printf("\n");


        printf("   --- Memory Information for device %d ---\n", i );
        //printf("Total global mem      :  %14.1f (bytes)\n", (double) prop.totalGlobalMem );
        //printf("Total global mem      :  %14.1f (kb)\n", prop.totalGlobalMem/1024.0);
        //printf("Total global mem      :  %14.1f (mb)\n", prop.totalGlobalMem/(1024.0*1024.0));
        printf("Total global mem      :  %14.1f (gb)\n",  prop.totalGlobalMem/(1024.0*1024.0*1024.0));
        printf("\n");

        printf( "   --- MP Information for device %d ---\n", i );
        printf( "Multiprocessor count :  %14d\n",
                    prop.multiProcessorCount );
        printf("Shared mem per mp     :  %14.1f (kb)\n", prop.sharedMemPerBlock/1024. );
        printf("Registers per mp      :  %14.1f (kb)\n", prop.regsPerBlock/1024. );
        printf("Threads in warp       :  %14d\n", prop.warpSize );
        printf("Max threads per block :  %14d\n",
                    prop.maxThreadsPerBlock );
        printf( "Max thread dimensions:  (%d, %d, %d)\n",
                    prop.maxThreadsDim[0], prop.maxThreadsDim[1],
                    prop.maxThreadsDim[2] );
        printf( "Max grid dimensions  :  %d, %d, %d\n",
                    prop.maxGridSize[0], prop.maxGridSize[1],
                    prop.maxGridSize[2] );
        printf( "\n" );
    }
}
