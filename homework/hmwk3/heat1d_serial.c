#include <stdio.h>
#include <stdlib.h>

#include <math.h>

/* Initial condition */
double init(double x)
{
    return exp(-5*x*x);;
}

int main(int argc, char** argv)
{
    /* ------------------------------ Input parameters -------------------------------- */
    int N = atoi(argv[1]);
    int nout = atoi(argv[2]);    

    /* Always print first and last time step, at a minimum */
    nout = (nout < 2) ? 2 : nout;  

    double L = 1;
    double Tfinal = 0.5;

    /* --------------------------- Numerical parameters ------------------------------- */
    double a = -L;
    double b = L;

    double dx = (b-a)/N;
    double dx2 = dx*dx;

    /* Write out meta data  */
    FILE *fout = fopen("heat1d.out","w");        
    fwrite(&N,1,sizeof(int),fout);
    fwrite(&a,1,sizeof(double),fout);
    fwrite(&b,1,sizeof(double),fout);

    /* ---------------------------- Initialize solution ------------------------------- */

    double *xmem = (double*) malloc((N+3)*sizeof(double));
    double *qmem = (double*) malloc((N+3)*sizeof(double));

    /* Offset indexing so we can index ghost cells with -1 */
    double *x = &xmem[1];
    double *q = &qmem[1];

    for(int i = -1; i <= N+1; i++)
    {
        x[i] = a + i*dx;
        q[i] = init(x[i]);
    }

    /* ----------------------------- Compute time step ---------------------------------*/
    /* Compute a stable time step
       1.  Estimate a stable time step 'dt_stable'.   This may not evenly divide Tfinal. 
       2.  Compute a minimum number M of time steps we need to take.
       3.  Divide Tfinal by M to get get a dt that is guaranteed smaller than dt_est and  
           satisfies M*dt = Tfinal.
    */
    double dt_stable = 0.95*dx2/2;        /* Stable time step */
    int M = ceil(Tfinal/dt_stable) + 1;   /* Compute M to guarantee we hit Tfinal */
    double dt = Tfinal/M;                 /* dt <= dt_stable; M*dt = Tfinal */

    /* More meta data */
    fwrite(&M,1,sizeof(int),fout);

    /* 
       Set up an array of 'n' values that tell us when to save our solution 
       so that we save exactly nout time steps at roughly equally spaced times.  
    */
    int *noutsteps = malloc(nout*sizeof(int));    
    double dM = ((double) M-1)/(nout-1);
    dM = (dM < 1) ? 1 : dM;
    for(int m = 0; m <= nout-1; m++)
    {
        noutsteps[m] = (int) floor(m*dM);
        //printf("%5d %5d\n",m,noutsteps[m]);
    }

    /* Output initial condition */
    double t = 0;
    int k = 0;
    //printf("Output frame %3d at t = %12.5e\n",k,t);
    fwrite(&t,1,sizeof(double),fout);       
    fwrite(&q[0],N+1,sizeof(double),fout);
    k++;  /* Number of output files created */

    /* --------------------------- Start time stepping ---------------------------------*/
    /* Store q^{n+1} */
    double *qpmem = (double*) malloc((N+3)*sizeof(double));
    double *qp = &qpmem[1];

    /* Time loop;  compute q^{n+1} at each time step  */
    for(int n = 0; n <= M-1; n++)
    {
        t += dt;

        /* No-flux boundary conditions */
        q[-1] = q[1];
        q[N+1] = q[N-1];
        for(int i = 0; i <= N; i++)
        {
            qp[i] = q[i] + dt*((q[i-1] - 2*q[i] + q[i+1])/dx2);            
        }
        if (n == noutsteps[k])
        {            
            //printf("Output frame %3d at t = %12.5e\n",k,t);
            fwrite(&t,1,sizeof(double),fout);       
            fwrite(&qp[0],N+1,sizeof(double),fout);
            k++;
        }

        q = qp;  
    }
    fclose(fout);

    free(xmem);
    free(qmem);
    free(qpmem);
    free(noutsteps);

    return 0;
}