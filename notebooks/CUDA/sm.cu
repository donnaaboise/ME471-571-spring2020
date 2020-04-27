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


int** allocate_2d(int N, int M)
{
    int rows = N;
    int cols = M;
    int   *qmem = (int*) malloc(rows*cols*sizeof(double));
    int **qrows = (int**) malloc(rows*sizeof(double*));

    for(int i = 0; i < rows; i++)
    {
        qrows[i] = &qmem[cols*i];
    }    
    return &qrows[0];
}

// Taken from "clockrate, computed in main

//#define CLOCK_RATE 745000      // On K40c (node 6 on Redhawk)
//#define CLOCK_RATE 875500      // On GTX Titan (node 5 on Redhawk)

#define CLOCK_RATE 1076000       // On Titan X (nodes 1-4 on Redhawk)

double cpuSecond()
{
    struct timeval tp;
    gettimeofday(&tp,NULL);
    return (double) tp.tv_sec + (double)tp.tv_usec*1e-6;
}

__device__ uint get_smid(void) {

     uint ret;
     asm("mov.u32 %0, %smid;" : "=r"(ret) );
     return ret;
}

__device__ void sleep(float t)
{    
    clock_t t0 = clock64();
    clock_t t1 = t0;
    while ((t1 - t0)/(CLOCK_RATE*1000.0f) < t)
    {
        t1 = clock64();
    }
}

__global__ void mykernel(float *t,uint *s) 
{
    int id = blockIdx.x;
    s[id] = (int) get_smid();
    sleep(t[id]);    
}

int main(int argc, char* argv[]) 
{
    cudaDeviceProp  prop;
    cudaGetDeviceProperties(&prop, 0); /* Only look at first GPU */

    int mp = prop.multiProcessorCount;
    clock_t clock_rate = prop.clockRate;
    printf("clock rate = %d\n",clock_rate);

    if (argc < 1)
    {
        printf("Usage : sm <num_blocks> <threads/blocks> <shared_memory> <scale_factor>\n");
        printf("    -- shared memory and scale factors are optional\n");
    }

    /* Read in number of blocks to launch */
    int num_blocks = atoi(argv[1]);    // number of blocks
    int threads_per_block = atoi(argv[2]);    // threads per block

    int shared_memory_per_thread = 0;
    if (argc > 3)
        shared_memory_per_thread = atoi(argv[3]);

    double scale_factor = 1;
    if (argc > 4)
        scale_factor = atof(argv[4]);

    int prt = 1;
    if (argc > 5)
        prt = atof(argv[5]);

    /* Create verctor of "work" for each block */
    float* t = (float*) malloc(num_blocks*sizeof(float));

    /* Create vector identifing which SM each block was loaded onto */
    uint* sm_id = (uint*) malloc(num_blocks*sizeof(uint));

    /* Allocate memory on the device */

    float *dev_t; 
    uint  *dev_sm_id;
    cudaMalloc( (void**) &dev_t,     num_blocks*sizeof(float));
    cudaMalloc( (void**) &dev_sm_id, num_blocks*sizeof(uint));

    printf("Memory requirement : %0.2f (kB)\n",num_blocks*(sizeof(float) + sizeof(uint))/(1024.0));

    /* thread work : Use a scaling factor to reduce the time; scale back once we are done */
    for(int i = 0; i < num_blocks; i++)
            t[i] = 1.0/scale_factor;            

    cudaMemcpy(dev_t, t, num_blocks*sizeof(float), cudaMemcpyHostToDevice);

    dim3 block(threads_per_block);
    dim3 grid(num_blocks);  /* N blocks */

    double start = cpuSecond();
    int shared_memory = shared_memory_per_thread*threads_per_block;

    mykernel<<<grid,block,shared_memory>>>(dev_t, dev_sm_id);
    CHECK(cudaDeviceSynchronize());
    double etime = cpuSecond() - start;

    CHECK(cudaPeekAtLastError());


    /* Copy contents of dev_t back to t */
    cudaMemcpy(sm_id, dev_sm_id, num_blocks*sizeof(uint), cudaMemcpyDeviceToHost);

    /* Post process data */
    int* blocks_per_SM = (int*) malloc(mp*sizeof(int));
    float* time_per_SM = (float*) malloc(mp*sizeof(float));
    printf("Device has %d SMs\n",mp);

    for(int i = 0; i < mp; i++)
    {
        blocks_per_SM[i] = 0;
        time_per_SM[i] = 0;
    }

    int **sm_list = (int**) allocate_2d(mp,num_blocks);

    printf("Distribution of blocks on SMs\n");
    for(int k = 0; k < num_blocks; k++)
    {
        int sm_num = sm_id[k];
        sm_list[sm_num][blocks_per_SM[sm_num]] = k;
        time_per_SM[sm_num] += scale_factor*t[k];
        blocks_per_SM[sm_num] += 1;
    }    


    printf("------------------------------------------------------------------------------\n");
    printf("\n");
    printf("SM    #blocks/SM  #threads/SM    #warps/SM     work/SM(s)       block list/SM\n");
    printf("------------------------------------------------------------------------------\n");
    for(int i = 0; i < mp; i++)
    {
        int bpSM = blocks_per_SM[i];
        int tpSM = blocks_per_SM[i]*threads_per_block;  
        int wpSM = blocks_per_SM[i]*(threads_per_block + 32-1)/32;
        printf("%2d %13d %12d %12d %10.2f  ",
               i,bpSM, tpSM, wpSM, time_per_SM[i]);

        if (prt)
        {
            if (bpSM > 0)
                printf("   (");

            for(int k = 0; k < bpSM; k++)
            {
                printf("%3d",sm_list[i][k]);
                if (k < bpSM-1)
                    printf(", ");
            }
            if (bpSM > 0)
                printf(")\n");
            else
                printf("\n");
        }
        else
        {
            printf("\n");
        }

    }
    printf("------------------------------------------------------------------------------\n");
    printf("\n");

    int total_threads = block.x*grid.x;
    printf("%32s %12d\n", "Threads per block",block.x*block.y);
    printf("%32s %12d\n", "Total number of blocks",grid.x);
    printf("%32s %12d\n", "Shared memory per block (bytes)",shared_memory);
    printf("%32s %12d\n", "Total number of threads",total_threads);
    printf("%32s %12g\n", "Time scaling factor",scale_factor);
    if (scale_factor != 1)
        printf("%27s %12.3f (s)\n","GPU Kernel Time (scaled)", scale_factor*etime);
    else
        printf("%27s %12.3f (s)\n","GPU Kernel Time", etime*scale_factor);
    printf("\n");

    cudaFree(dev_t);
    cudaFree(dev_sm_id);

    free(t);
    free(sm_id);
    free(blocks_per_SM);

    cudaDeviceReset();
}


