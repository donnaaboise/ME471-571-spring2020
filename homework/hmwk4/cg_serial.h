#ifndef CG_SERIAL_H
#define CG_SERIAL_H


#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif


typedef void (*matmult_t)(int N, double* x,double *Ax);

void cg_1d(int N, int mbc, int kmax, double tol, int prt, matmult_t matmult, 
        double *x, double *F, int* it_cnt);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif


#endif