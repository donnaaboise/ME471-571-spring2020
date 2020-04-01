#include <stdio.h>
#include <demo_util.h>

#define CLOCK_RATE 1124000     // in kHZ

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

__global__ void worker(float *t,int *s) 
{
    int id = blockIdx.x;
    s[id] = get_smid();
    sleep(t[id]);
}

#define N 1536

int main(void) 
{
    float *dev_t;
    int *dev_s;
    float t[N];
    int s[N];
    float SM[8] = {0,0,0,0,0,0,0,0};
    int i;

    /* Allocate memory on the device */
    cudaMalloc( (void**)&dev_t, N*sizeof(float));
    cudaMalloc( (void**)&dev_s, N*sizeof(int));

    random_seed();

    float maxt = 0;
    for(i = 0; i < N; i++)
    {           
        t[i] = 1; // 10*random_number();
        maxt = t[i] > maxt ? t[i] : maxt;
    }

    cudaMemcpy(dev_t, t, N*sizeof(float), cudaMemcpyHostToDevice);

    worker<<<48,32>>>(dev_t, dev_s);

    /* Copy contents of dev_t back to t */
    cudaMemcpy(s, dev_s, N*sizeof(float), cudaMemcpyDeviceToHost);

    for(i = 0; i < N; i++)
    {
        printf( "Block %2d worked for %8.4f seconds on SM %d\n",i,t[i],s[i]);
        SM[s[i]] += t[i];
    }
    printf("Max t = %8.3f\n",maxt);
    printf("\n");
    for(i = 7; i>= 0; i--)
    {
        printf("SM[%d] = %8.4f\n",i,SM[i]);
    }
    cudaFree(dev_t);

}


