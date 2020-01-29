#ifndef DEMO_UTIL_H
#define DEMO_UTIL_H

#include <time.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

typedef enum
{
    SILENT = 0,
    ESSENTIAL,
    PRODUCTION,
    INFO,
    DEBUG
} loglevel_t;


/* Array handling routines.  These allocate or delete memory */
void empty_array(int n,double **x);
void zeros_array(int n,double **x);
void ones_array(int n,double **x);
void constant_array(int n,double **x, double value);
void random_array(int n, double **array);
void linspace_array(double a,double b,int n,double **x);
void delete_array(void **x);

void empty_array2(int nrows,int ncols, double ***x);
void delete_array2(void ***A);

void char_array(int n, char **c);
void delete_char_array(char **c);

/* Operations on arrays */
double sum_array(int n, double *x);

/* I/O routines */
void read_int   (int argc, char** argv, const char arg[], int* value, int *err);
void read_double(int argc, char** argv, const char arg[], double* value,int *err);
void read_string(int argc, char** argv, const char arg[], char* value,int *err);


void read_loglevel(int argc, char** argv);
void print_global(const char* format, ... );
void print_essential(const char* format, ... );
void print_production(const char* format, ... );
void print_info(const char* format, ... );
void print_debug(const char* format, ... );

/* Miscellaneous */
void set_rank(int my_rank);
double random_number();
void random_seed();
int pow2(int p);
int random_int(int m, int n);
void sleep(double t_total);



#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
