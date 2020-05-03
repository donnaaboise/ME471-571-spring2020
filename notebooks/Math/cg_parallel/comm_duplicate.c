#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

#include "solver.h"

void get_endpoints_indices(int N_global, int *i1, int* i2, int *j1, int*j2, 
                           MPI_Comm comm_cart)
{
    int maxdims = 2;
    int dims[2],periods[2], mycoords[2];
    MPI_Cart_get(comm_cart,maxdims,dims,periods,mycoords);

    int Nx = N_global/dims[0];
    int Ny = N_global/dims[1];

    bctype_t bc_type[4];
    get_phys_bdry(bc_type, comm_cart);

    *i1 = (bc_type[0] == INTERIOR) ? 0  : 1;
    *j1 = (bc_type[2] == INTERIOR) ? 0  : 1;
    *i2 = (bc_type[1] == INTERIOR) ? Nx : Nx-1;
    *j2 = (bc_type[3] == INTERIOR) ? Ny : Ny-1;
}


void get_phys_bdry(bctype_t bc_type[], MPI_Comm comm_cart)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int maxdims = 2;
    int dims[2],periods[2], mycoords[2];
    MPI_Cart_get(comm_cart,maxdims,dims,periods,mycoords);

    int source, dest, dir, disp, ierr;

    /* left edge */
    dir = 0;
    disp = 1;
    ierr = MPI_Cart_shift(comm_cart, dir, disp, &source, &dest);
    bc_type[0] = (source == MPI_PROC_NULL) ? PHYS : INTERIOR;

    /* right edge */
    dir = 0;
    disp = -1;
    ierr = MPI_Cart_shift(comm_cart, dir, disp, &source, &dest);
    bc_type[1] = (source == MPI_PROC_NULL) ? PHYS : INTERIOR;

    /* bottom edge */
    dir = 1;
    disp = 1;
    ierr = MPI_Cart_shift(comm_cart, dir, disp, &source, &dest);
    bc_type[2] = (source == MPI_PROC_NULL) ? PHYS : INTERIOR;

    /* top edge */
    dir = 1;
    disp = -1;
    ierr = MPI_Cart_shift(comm_cart, dir, disp, &source, &dest);
    bc_type[3] = (source == MPI_PROC_NULL) ? PHYS : INTERIOR;
}


void comm(int N_global, double **u, MPI_Comm comm_cart)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int maxdims = 2;
    int dims[2],periods[2], mycoords[2];
    MPI_Cart_get(comm_cart,maxdims,dims,periods,mycoords);

    int Nx = N_global/dims[0];
    int Ny = N_global/dims[1];

    double *ghost_buff[2];
    ghost_buff[0] = malloc((N_global+1)*sizeof(double));
    ghost_buff[1] = malloc((N_global+1)*sizeof(double));

    /* Get indices endpoints for unknowns on this patch */
    int i1,i2,j1,j2;
    get_endpoints_indices(N_global, &i1, &i2, &j1, &j2,comm_cart);

    bctype_t bc_type[4];
    get_phys_bdry(bc_type, comm_cart);
    int tag = 0;

    int source, dest;

    /* Shift right (fill left ghost cells) : Send Nx-1 from left grid to -1
       of the right grid */
    {        
        int disp = 1;
        int dir = 0;
        for(int j = 0; j < Ny+1; j++)
            ghost_buff[0][j] = u[Nx-1][j];

        int ierr = MPI_Cart_shift(comm_cart, dir, disp, &source, &dest);

        MPI_Sendrecv(ghost_buff[0], Ny+1, MPI_DOUBLE, dest, tag, 
                     ghost_buff[1], Ny+1, MPI_DOUBLE, source, tag, 
                     comm_cart, MPI_STATUS_IGNORE);

        if (bc_type[0] == INTERIOR)
            for(int j = 0; j < Ny+1; j++)
                u[-1][j]  = ghost_buff[1][j];       
    }     

    /* Shift left (fill right ghost cells) : Send 1 from right grid to Nx+1
       of the left grid */
    {        
        int disp = -1;
        int dir = 0;
        for(int j = 0; j < Ny+1; j++)
            ghost_buff[0][j] = u[1][j];

        int ierr = MPI_Cart_shift(comm_cart, dir, disp, &source, &dest);

        MPI_Sendrecv(ghost_buff[0], Ny+1, MPI_DOUBLE, dest, tag, 
                     ghost_buff[1], Ny+1, MPI_DOUBLE, source, tag, 
                     comm_cart, MPI_STATUS_IGNORE);


        if (bc_type[1] == INTERIOR)
            for(int j = 0; j < Ny+1; j++)
                u[Nx+1][j]  = ghost_buff[1][j];

    }

    /* Shift up (fill bottom ghost cells) : Send Ny-1 from bottom grid to -1
       of the top grid */
    {
        int disp = 1;
        int dir = 1;
        for(int i = 0; i < Nx+1; i++)
            ghost_buff[0][i] = u[i][Ny-1];

        int ierr = MPI_Cart_shift(comm_cart, dir, disp, &source, &dest);

        MPI_Sendrecv(ghost_buff[0], Nx+1, MPI_DOUBLE, dest, tag, 
                     ghost_buff[1], Nx+1, MPI_DOUBLE, source, tag, 
                     comm_cart, MPI_STATUS_IGNORE);

        if (bc_type[2] == INTERIOR)
            for(int i = 0; i < Nx+1; i++)
                u[i][-1]  = ghost_buff[1][i];
    }

    /* Shift down (fill top ghost cells) : Send 1 from top grid to Ny+1
       of the bottom grid */
    {        
        int disp = -1;
        int dir = 1;
        for(int i = 0; i < Nx+1; i++)
            ghost_buff[0][i] = u[i][1];

        int ierr = MPI_Cart_shift(comm_cart, dir, disp, &source, &dest);

        MPI_Sendrecv(ghost_buff[0], Nx+1, MPI_DOUBLE, dest, tag, 
                     ghost_buff[1], Nx+1, MPI_DOUBLE, source, tag, 
                     comm_cart, MPI_STATUS_IGNORE);

        if (bc_type[3] == INTERIOR)
            for(int i = 0; i < Nx+1; i++)
                u[i][Ny+1]  = ghost_buff[1][i];                

    }
    free(ghost_buff[0]);
    free(ghost_buff[1]);
}