#include <mpi.h>
#include <math.h>
#include <stdio.h>  
#include <stddef.h>  /* offsetof */
#include <stdlib.h>

#define PI 3.14159265358979323846264338327

void mpi_debug();

double utrue(double x, double y)
{
    double u;
    double pi2;
    pi2 = 2*PI;
    //u = cos(pi2*x)*sin(PI*y);

    double r2 = pow((x - 0.15),2) + pow((y-0.65),2);
    u = exp(-50*r2);
    return u;
}


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

typedef struct
{
    double t; 
    MPI_Datatype localarray_t;
} struct_timeinfo_t;

void create_timeinfo_type(MPI_Datatype localarray_t, MPI_Datatype *timeinfo_t)
{
    int blocksize = 2;
    int block_lengths[2] = {1,1}; /* Use blocksize to dimension array */


    /* Set up types */
    MPI_Datatype typelist[2];
    typelist[0] = MPI_DOUBLE;
    typelist[1] = localarray_t;

    /* Set up displacements */
    MPI_Aint disp[4];
    disp[0] = offsetof(struct_timeinfo_t,t);
    disp[1] = offsetof(struct_timeinfo_t,localarray_t);

    MPI_Type_create_struct(blocksize,block_lengths, disp, typelist, timeinfo_t);
    MPI_Type_commit(timeinfo_t);
}



typedef struct 
{
    double xgrid[2];
    double ygrid[2];
    int N[2];
    int Mout;
} struct_header_t;

void create_header_type(MPI_Datatype *header_t)
{
    int blocksize = 4;
    int block_lengths[4] = {2,2,2,1}; /* Use blocksize to dimension array */

    /* Set up types */
    MPI_Datatype typelist[4];
    typelist[0] = MPI_DOUBLE;
    typelist[1] = MPI_DOUBLE;
    typelist[2] = MPI_INT;
    typelist[3] = MPI_INT;

    /* Set up displacements */
    MPI_Aint disp[4];
    disp[0] = offsetof(struct_header_t,xgrid);
    disp[1] = offsetof(struct_header_t,ygrid);
    disp[2] = offsetof(struct_header_t,N);
    disp[3] = offsetof(struct_header_t,Mout);

    MPI_Type_create_struct(blocksize,block_lengths, disp, typelist, header_t);
    MPI_Type_commit(header_t);
}


void main(int argc, char** argv)
{
   /* ------------------------------- MPI Initialization ------------------------------ */
    MPI_Init(&argc, &argv);
 
 
    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);



    /* ------------------------------- Read command line ------------------------------ */
    int p = atoi(argv[1]);
    int prows = atoi(argv[2]);   /* Number of procs in x direction */
    int pcols = atoi(argv[3]);   /* Number of procs in y direction */
    int Mout = atoi(argv[4]);
    // int debug = atoi(argv[5]);

#if 0
    if (debug)
    {
        mpi_debug();        
    }
#endif    

    if (prows*pcols != nprocs)
    {
        printf("prows*pcols != nprocs\n");
        exit(0);
    }

    /* ------------------------------- Other user parameters -------------------------- */
    /* header is [ax,bx] x [ay,by] */
    double ax = 0;
    double bx = 1;
    double ay = 0;
    double by = 1;

    /* Number of mesh cells in each direction */
    int Nx = pow(2,p);
    int Ny = Nx;

    /* --------------------------------- Communicator --------------------------------- */

    MPI_Comm comm_cart;
    int ndim = 2;
    int dims[2] = {prows,pcols};  /* (0,0) is lower left corner */
    int periodicity[2] = {0,0};
    int reorder = 0;
    MPI_Cart_create(MPI_COMM_WORLD, ndim, dims,  periodicity, reorder, &comm_cart);


    int mycoords[2];
    MPI_Cart_coords(comm_cart,rank,2,mycoords);

    /* ------------------------------- Numerical parameters --------------------------- */

    int Nx_local = Nx/dims[0];
    int Ny_local = Ny/dims[1];

    double w[2] = {(bx-ax)/dims[0], (by-ay)/dims[1]};

    double dx = (bx-ax)/Nx;
    double dy = (by-ay)/Ny;

    /* --------------------------------- Header ----------------------------------------*/
    MPI_File file;
    MPI_File_open(MPI_COMM_WORLD, "gridsoln.out", 
                  MPI_MODE_CREATE|MPI_MODE_WRONLY,
                  MPI_INFO_NULL, &file);

    MPI_Datatype header_t;
    create_header_type(&header_t);

    struct_header_t header;    
    if (rank == 0)
    {
        header.xgrid[0] = ax;
        header.xgrid[1] = bx;
        header.ygrid[0] = ay;
        header.ygrid[1] = by;
        header.N[0] = Nx;
        header.N[1] = Ny;
        header.Mout = Mout;

        MPI_File_write(file,&header,1,header_t, MPI_STATUS_IGNORE);  
    }

    /* --------------------- Create view for this processor into file ------------------*/

    /* Create a file type : Data is distributed in a grid across processors */
    int ndims = 2;
    int globalsize[2] = {Nx+1,Ny+1};
    int localsize[2] = {Nx_local+1,Ny_local+1};
    int starts[2] = {mycoords[0]*Nx_local, mycoords[1]*Ny_local};
    int order = MPI_ORDER_C;  /* MPI_ORDER_FORTRAN */

    MPI_Datatype localarray_t;
    MPI_Type_create_subarray(ndims, globalsize, localsize, starts, order, 
                             MPI_DOUBLE, &localarray_t);

    MPI_Type_commit(&localarray_t);


    MPI_Datatype timeinfo_t;
    create_timeinfo_type(localarray_t,&timeinfo_t);

    MPI_Aint extent;
    MPI_Type_extent(header_t,&extent); 
    MPI_Offset disp = extent;   /* Should be 48 bytes (4*8 + 3*4 + 4(extra)) */
    MPI_File_set_view(file, disp,  MPI_DOUBLE, localarray_t, 
                           "native", MPI_INFO_NULL);


    /* ----------------------------------- Initialize data ---------------------------- */


    double **const u = allocate_2d(Nx_local,Ny_local,0);

    for(int j = 0; j <= Ny_local; j++)
    {
        double y = ay + w[1]*mycoords[1] + j*dy;
        for(int i = 0; i <= Nx_local; i++)
        {
            double x = ax + w[0]*mycoords[0] + i*dx;
            u[i][j] = utrue(x,y);            
        }
    }

    /* ------------------------------------ "Time step" ------------------------------- */

    /* Write out initial time and solution */
    double t = 0;

    MPI_Offset offset;
    MPI_Aint extent_la;
    MPI_Type_extent(localarray_t,&extent_la); 

    /* Write out the time and reposition the file handle */
    if (rank == 0)
    {
        MPI_File_write(file,&t, 1, MPI_DOUBLE,MPI_STATUS_IGNORE);                        
    }
    disp += sizeof(double);
    MPI_File_set_view(file, disp,  MPI_DOUBLE, localarray_t, 
                      "native", MPI_INFO_NULL);

    /* Output local array and update displacement */
    int usize = (Nx_local+1)*(Ny_local+1);
    MPI_File_write_all(file, &u[0][0], usize, MPI_DOUBLE, MPI_STATUS_IGNORE);
    disp += extent_la;

    double dt = 0.1;
    for(int n = 1; n < Mout; n++)
    {
        t += dt;
        for(int j = 0; j <= Ny_local; j++)
        {
            for(int i = 0; i <= Nx_local; i++)
            {
                u[i][j] *= exp(-dt);   /* Fake heat diffusion */
            }
        }
        /* Write out the time and reposition the file handle */
        if (rank == 0)
        {
            MPI_File_write(file,&t, 1, MPI_DOUBLE,MPI_STATUS_IGNORE);                        
        }
        disp += sizeof(double);
        MPI_File_set_view(file, disp,  MPI_DOUBLE, localarray_t, 
                          "native", MPI_INFO_NULL);

        /* Collective write; update file displacement */
        MPI_File_write_all(file, &u[0][0], usize, MPI_DOUBLE, MPI_STATUS_IGNORE);
        disp += extent_la;
    }


    /* ------------------------------- Clean up --------------------------------------- */

    MPI_Type_free(&localarray_t);
    MPI_Type_free(&header_t);

    delete_2d(0,(double***) &u);

    MPI_File_close(&file);

    MPI_Finalize();

}