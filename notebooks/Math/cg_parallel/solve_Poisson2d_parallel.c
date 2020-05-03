#include <math.h>
#include <stdlib.h>   /* For atoi */
#include <stdio.h>
#include <string.h>
#include <stddef.h>   /* offsetof */

#include <mpi.h>

#include "solver.h"   /* For f_t def */

#define PI 3.14159265358979323846264338327

typedef double (*f_t)(double x, double y);


double** allocate_2d(int N, int M, int mbc)
{
    int rows = N + 1 + 2*mbc;
    int cols = M + 1 + 2*mbc; 

    double   *qmem = malloc(rows*cols*sizeof(double));
    double **qrows = malloc(rows*sizeof(double*));

    for(int i = 0; i < rows; i++)
    {
        qrows[i] = &qmem[cols*i + mbc];
    }    
    return &qrows[mbc];
}

void delete_2d(int mbc, double ***q)
{
    free(&(*q)[-mbc][-mbc]);
    free(&(*q)[-mbc]);
    *q = NULL;
}


double utrue0(double x, double y)
{
    return x*x*x;
}


double rhs0(double x,double y)
{
    return 6*x;
}


double utrue1(double x,double y)
{
    double pi2 = 2*PI;
    double u = sin(pi2*x)*sin(pi2*y);
    return u;
}


double rhs1(double x,double y)
{
    double pi2 = 2*PI;
    double lapf = -2*(pi2)*(pi2)*utrue1(x,y);
    return lapf;
}


double utrue2(double x, double y)
{
    double pi = PI;
    double a = 20.0;
    double pi2 = 2*PI;

    double x1 = x-0.35;
    double y1 = y-0.55;
    double r2 = x1*x1 + y1*y1;
    double f = exp(-a*r2);
    double g = sin(pi2*x)*sin(pi2*y);

    double u = f*g;
    return u;
}

static
double dot2(double u[2], double v[2])
{
    return u[0]*v[0] + u[1]*v[1];
}

double rhs2(double x, double y)
{
    double pi = PI;
    double a = 20;
    double pi2 = 2*PI;

    double x1 = x-0.35;
    double y1 = y-0.55;
    double r2 = x1*x1 + y1*y1;

    double f = exp(-a*r2);
    double g = sin(pi2*x)*sin(pi2*y);

    double fx  = -2*a*x1*f;
    double fxx = -2*a*(x1*fx + f);
    double fy  = -2*a*y1*f;
    double fyy = -2*a*(y1*fy + f);

    double lapf = fxx + fyy;
    double gradf[2] = {fx, fy};


    double gx  = pi2*cos(pi2*x)*sin(pi2*y);
    double gxx = -pow(pi2,2.0)*sin(pi2*x)*sin(pi2*y);
    double gy  = pi2*sin(pi2*x)*cos(pi2*y);
    double gyy = -pow(pi2,2.0)*sin(pi2*x)*sin(pi2*y);

    double lapg = gxx + gyy;
    double gradg[2] = {gx, gy};

    double u = f*g;
    double lapu = g*lapf + 2*dot2(gradg,gradf) + f*lapg;
    return lapu;
}


typedef struct 
{
    double xgrid[2];
    double ygrid[2];
    int N[2];
} struct_header_t;

void create_header_type(MPI_Datatype *header_t)
{
    int blocksize = 3;
    int block_lengths[3] = {2,2,2}; /* Use blocksize to dimension array */

    /* Set up types */
    MPI_Datatype typelist[3];
    typelist[0] = MPI_DOUBLE;
    typelist[1] = MPI_DOUBLE;
    typelist[2] = MPI_INT;

    /* Set up displacements */
    MPI_Aint disp[3];
    disp[0] = offsetof(struct_header_t,xgrid);
    disp[1] = offsetof(struct_header_t,ygrid);
    disp[2] = offsetof(struct_header_t,N);

    MPI_Type_create_struct(blocksize,block_lengths, disp, typelist, header_t);
    MPI_Type_commit(header_t);
}


static double s_dx;
static double s_dy;

static int s_inh_flag;

void matmult(int N_global, double** x, double **Ax, MPI_Comm comm_cart)
{
    /* communicate ghost cells */
    comm(N_global,x,comm_cart);

    /* Get indices endpoints for unknowns on this patch */
    int i1,i2,j1,j2;
    get_endpoints_indices(N_global, &i1, &i2, &j1, &j2,comm_cart);

    double dx2 = s_dx*s_dx;
    double dy2 = s_dx*s_dx;

    /* Compute A*x */
    for(int i = i1; i <= i2; i++)
        for(int j = j1; j <= j2; j++)
        {
            double uxx = (x[i-1][j] - 2*x[i][j] + x[i+1][j])/dx2;
            double uyy = (x[i][j-1] - 2*x[i][j] + x[i][j+1])/dy2;
            Ax[i][j] = uxx + uyy;                      
        }
}


int main(int argc, char** argv)
{

    MPI_Init(&argc, &argv);
    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* ----------------------------------- user defined parameters -------------------- */
    int        dimx = atoi(argv[1]);
    int        dimy = atoi(argv[2]);
    int    N_global = atoi(argv[3]);
    double      tol = atof(argv[4]);    // 1e-12
    int        kmax = atoi(argv[5]);   // 1000
    int         prt = atoi(argv[6]);

#if 0    
    int  rhs_choice = atoi(argv[7]);
    method_t mth    = atoi(argv[8]);
    int debug       = atoi(argv[9]);
#else    
    /* Allow optional arguments */
    int rhs_choice;    /* Choose different RHS to see diferent convergence behaviors */
    int debug; 
    method_t method;
    if (argc < 10)
    {
        debug = 0;     /* No debugging by default */
        if (argc < 9)
        {
            method = CG;        /* Use an 'enum' here - see cg_parallel.h */
            if (argc < 8)
                rhs_choice = 1;    /* Use RHS from homework problem */
            else
                rhs_choice = atoi(argv[7]);
        }        
        else
        {
            rhs_choice = atoi(argv[7]);
            method = atoi(argv[8]);
        }
    }
    else
    {
        rhs_choice = atoi(argv[7]);
        method = atoi(argv[8]);
        debug = atoi(argv[9]);        
    }

#endif    

    /* ------------------------- Post process input arguments ------------------------- */

    f_t utrue, rhs;
    switch(rhs_choice)
    {
        case 1:
            /* Use RHS from homework problem */
            utrue = utrue1;  
            rhs = rhs1;
            break;
        case 2:
            /* Use RHS with exp. and trig functions */
            utrue = utrue2;  
            rhs = rhs2;
            break;
        default:
            printf("Error : No valid RHS choice specified\n");
            exit(0);
    }

    int dims[2] = {dimx,dimy};


    /* ----------------------------- Print out user options --------------------------- */
    if (rank == 0)
    {
        printf("\n");
        printf("Input arguments\n");
        printf("%15s %10d\n","Dim[0]",dims[0]);
        printf("%15s %10d\n","Dim[1]",dims[1]);
        printf("%15s %10d\n","N",N_global);
        printf("%15s %10.0e\n","tol",tol);
        printf("%15s %10d\n","kmax",kmax);
        printf("%15s %10d\n","prt",prt);
        printf("\n");
        printf("Optional arguments\n");
        printf("%15s %10d   (1=utrue1; 2=utrue2)\n","Problem",rhs_choice);
        printf("%15s %10d   (CG=2; BiCG=3)\n","Method",method);
        printf("%15s %10d\n","debug",debug);        
        printf("\n");
        printf("\n");
        printf("Iterations\n");
    }



    /* ------------------------------- Start debugging -------------------------------- */

    if (debug == 1)
        mpi_debug();        

    /* ----------------------------- Create Cartesian Topology ------------------------ */

    /* Fill in q[-1] and q[1] with data from neighbors */
    MPI_Comm comm_cart;
    int ndim = 2;
    int periodicity[2] = {0,0};
    int reorder = 0;
    MPI_Cart_create(MPI_COMM_WORLD, ndim, dims,  periodicity, reorder, &comm_cart);

    int mycoords[2];
    MPI_Cart_coords(comm_cart,rank,2,mycoords);

    bctype_t bc_type[4];
    get_phys_bdry(bc_type, comm_cart);

    int i1c, i2c, j1c, j2c;
    get_endpoints_indices(N_global, &i1c, &i2c, &j1c, &j2c,comm_cart);

    /* To be sure these don't get changed */
    const int i1 = i1c;
    const int i2 = i2c;
    const int j1 = j1c;
    const int j2 = j2c;

    /* -------------------------------- User parameters ------------------------------- */
    /* Number of ghost cells */
    int mbc = 1;  

    double a = 0;
    double b = 1;

    /* ----------------------------------- Numerical parameters ----------------------- */
    int Nx = N_global/dimx;
    int Ny = N_global/dimy;

    double dx = (b-a)/N_global;
    double dy = (b-a)/N_global;
    double w[2] = {(b-a)/dims[0], (b-a)/dims[1]};

    /* Set static variable to be used in matmult */
    s_dx = dx;   
    s_dy = dy;  

    /* --------------------------------- Header ----------------------------------------*/
    MPI_File file;
    MPI_File_open(MPI_COMM_WORLD, "iterative.out", 
                  MPI_MODE_CREATE|MPI_MODE_WRONLY,
                  MPI_INFO_NULL, &file);

    MPI_Datatype header_t;
    create_header_type(&header_t);

    struct_header_t header;    
    if (rank == 0)
    {
        header.xgrid[0] = a;
        header.xgrid[1] = b;
        header.ygrid[0] = a;
        header.ygrid[1] = b;
        header.N[0] = N_global;
        header.N[1] = N_global;

        MPI_File_write(file,&header,1,header_t, MPI_STATUS_IGNORE);  
    }

    /* --------------------- Create view for this processor into file ------------------*/

    /* Create a file type : Data is distributed in a grid across processors */
    int ndims = 2;
    int globalsize[2] = {N_global + 1,N_global + 1};
    int localsize[2] = {Nx+1,Ny+1};
    int starts[2] = {mycoords[0]*Nx, mycoords[1]*Ny};
    int order = MPI_ORDER_C;  /* MPI_ORDER_FORTRAN */

    MPI_Datatype localarray_t;
    MPI_Type_create_subarray(ndims, globalsize, localsize, starts, order, 
                             MPI_DOUBLE, &localarray_t);

    MPI_Type_commit(&localarray_t);

    MPI_Aint extent;
    MPI_Type_extent(header_t,&extent); 
    MPI_Offset disp = extent;   /* Should be 48 bytes (4*8 + 3*4 + 4(padding)) */
    MPI_File_set_view(file, disp,  MPI_DOUBLE, localarray_t, 
                           "native", MPI_INFO_NULL);

    /* -------------------------------------------- RHS ------------------------------- */

    int nsize = (Nx + 1 + 2*mbc)*(Ny + 1 + 2*mbc);

    double **const F = allocate_2d(Nx,Ny,mbc);
    memset(&F[-mbc][-mbc],0,nsize*sizeof(double));

    for(int i = i1; i <= i2; i++)
    {
        for(int j = j1; j <=j2; j++)
        {
            double x = a + w[0]*mycoords[0] + i*dx;
            double y = a + w[1]*mycoords[1] + j*dy;
            F[i][j] = rhs(x,y);
        }
    }
    double dx2 = dx*dx;
    double dy2 = dy*dy;

    /* Include non-zero boundary conditions */
    for(int j = j1; j <= j2; j++)
    {
        double y = a + w[1]*mycoords[1] + j*dy;
        if (bc_type[0] == PHYS)
            F[1][j] -= utrue(a,y)/dx2;

        if (bc_type[1] == PHYS)
            F[Nx-1][j] -= utrue(b,y)/dx2;
    }
    for(int i = i1; i <= i2; i++)
    {
        double x = a + w[0]*mycoords[0] + i*dx;
        if (bc_type[2] == PHYS)
            F[i][1]   -= utrue(x,a)/dy2;

        if (bc_type[3] == PHYS)
            F[i][Ny-1] -= utrue(x,b)/dy2;
    }

    /* -----------------------------------------------------------------------------------
      Solve Poisson Problem using CG, BICG or stationery methods (Jacobi, Gauss-Seidel)
    ----------------------------------------------------------------------------------- */

    double **const u = allocate_2d(Nx,Ny,mbc);
    memset(&u[-mbc][-mbc],0,nsize*sizeof(double));

    int it_cnt;
    double res;

    switch(method)
    {
        case CG:
            cg_2d(N_global,mbc,kmax, tol, prt, matmult,
                  &F[0],&u[0],&it_cnt, &res, comm_cart);
            break;
        case BICG:
            bicg_2d(N_global,mbc,kmax, tol, prt, matmult,
                  &F[0],&u[0],&it_cnt, &res, comm_cart);
            break;
        case JA:
        case GS:
            splitting(N_global,mbc,kmax, tol, prt, matmult,method, 
                      &F[0],&u[0],&it_cnt,&res, comm_cart);
            break;        
        default:        
            if (rank == 0)
                printf("No valid method\n");
    }

    /* ----------------------------------------------------------------
       Calculate error and report results
    ---------------------------------------------------------------- */

    /* Include boundary conditions in computed solution */
    for(int j = 0; j < Ny+1; j++)
    {
        double y = a + w[1]*mycoords[1] + j*dy;
        if (bc_type[0] == PHYS)
            u[0][j] = utrue(a,y);
        if (bc_type[1] == PHYS)
            u[Nx][j] = utrue(b,y);
    }

    for(int i = 0; i < Nx+1; i++)
    {
        double x = a + w[0]*mycoords[0] + i*dx;
        if (bc_type[2] == PHYS)
            u[i][0] = utrue(x,a);
        if (bc_type[3] == PHYS)
            u[i][Ny] = utrue(x,b);
    }

    /* Collective write: first copy data to array without ghost cells */
    int mbcout = 0;
    int usize = (Nx  + 1)*(Ny + 1);    
    double **const uout = allocate_2d(Nx,Ny,mbcout);    

    double err_local[3] = {0,0,0};
    for(int i = 0; i < Nx+1; i++)
        for(int j = 0; j < Ny+1; j++)
        {
            double x = a + w[0]*mycoords[0] + i*dx;
            double y = a + w[1]*mycoords[1] + j*dy;
            double udiff = u[i][j] - utrue(x,y);
            err_local[0] += fabs(udiff)*dx*dy;
            err_local[1] += fabs(udiff*udiff)*dx*dy;
            err_local[2] = fabs(udiff) > err_local[2] ? fabs(udiff) : err_local[2];    
            uout[i][j] = u[i][j];        
        }

    double err[3];
    MPI_Allreduce(err_local,err,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);        
    MPI_Allreduce(&err_local[2],&err[2],1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);        
    err[1] = sqrt(err[1]);    /* 2-norm */

    if (rank == 0)
    {
        printf("\n");
        printf("Results\n");
        printf("%10d %10d %16.8e %12.4e %12.4e %12.4e\n",N_global,it_cnt, res,
               err[0],err[1],err[2]);

        FILE *fout = fopen("error.dat","w");    
        fwrite(&N_global,1,sizeof(int),fout);
        fwrite(&it_cnt,1,sizeof(int),fout);
        fwrite(&res,1,sizeof(double),fout);
        fwrite(err,3,sizeof(double),fout);    
        fclose(fout);
    }

    /* ----------------------------------------------------------------
       Write out solution
    ---------------------------------------------------------------- */

    MPI_File_write_all(file, &uout[0][0], usize, MPI_DOUBLE, MPI_STATUS_IGNORE);

    MPI_Type_free(&localarray_t);
    MPI_Type_free(&header_t);
    MPI_File_close(&file);

    delete_2d(mbc,(double***) &F);
    delete_2d(mbc,(double***) &u);
    delete_2d(mbcout,(double***) &uout);

    MPI_Finalize();

    return 0;
}

