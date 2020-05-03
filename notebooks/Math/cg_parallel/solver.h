#ifndef CG_SERIAL_H
#define CG_SERIAL_H


#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

typedef enum method
{
    JA = 0,
    GS = 1,
    CG = 2,
    BICG = 3
} method_t;

typedef enum bctype
{
    PHYS = 0,
    INTERIOR = 1
} bctype_t;



double** allocate_2d(int N, int M, int mbc);
void delete_2d(int mbc, double ***q);

typedef void (*matmult2d_t)(int N, double** x, double **Ax, MPI_Comm comm_cart);

void comm(int N, double **u, MPI_Comm comm_cart);

void matmult(int N, double** x, double **Ax,MPI_Comm comm_cart);

void get_phys_bdry(bctype_t bc_type[], MPI_Comm comm_cart);

void get_endpoints_indices(int N_global, int *i1, int* i2, int *j1, int*j2, MPI_Comm comm_cart);


void cg_2d(int N, int mbc, int kmax, double tol, int prt, matmult2d_t matmult, 
        double **x, double **F, int* it_cnt, double *res, MPI_Comm comm_cart);

void bicg_2d(int N, int mbc, int kmax, double tol, int prt, matmult2d_t matmult, 
             double **x, double **F, int* it_cnt, double *res, MPI_Comm comm_cart);

void splitting(int N, int mbc, int kmax, double tol, int prt,  
               matmult2d_t matmult, method_t method, double **F, 
               double **u, int* it_cnt,double *res, MPI_Comm comm_cart);

void mpi_debug();

#ifdef __cplusplus
#if 0
{
#endif
}
#endif


#endif