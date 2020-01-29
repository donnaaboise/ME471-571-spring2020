#include "array_util.h"

#include <stdio.h>
#include <stdlib.h>

#include <math.h>
#include <time.h>

#include <stdarg.h>   /* for va_list, va_start, va_end, vprintf */
#include <string.h>   /* For atoi */

static int s_rank;
static loglevel_t s_loglevel;


/* --------------------------------------------------------------------
    Arrays : 

    Routines that allocate memory and possibly assign
    values.  Call these as follows : 

        double *x;    // No memory allocated yet
        int n = 10;

        empty_array(n,&x);  // x is now ready for use 
        x[0] = 1;
        x[1] = 2;
        // ...

        delete_array(&x);  // clean up when you are done!

    Array routine : 
        empty_array     : Allocate memory, but don't assign any values.
        ones_array      : Create an array of ones
        linspace_array  : Analogous to Matlab's 'linspace' 
        random_array    : Array of random numbers in [0,1]
        delete_array    : Deletes memory for x
   ------------------------------------------------------------------ */ 

static
void pointer_array(int n, void*** x)
{
    *x = malloc(n*sizeof(void*));
    if (x == NULL)
    {
        print_essential("pointer_array : array allocation error\n");
        exit(0);
    }
}

void empty_array(int n,double **x)
{
    *x = malloc(n*sizeof(double));
    if (x == NULL)
    {
        print_essential("empty_array : array allocation error\n");
        exit(0);
    }
}

void empty_array2(int nrows,int ncols, double ***A)
{
    double *x;
    int i;

    empty_array(nrows*ncols,&x);     
    pointer_array(nrows,(void***) A); 
    for(i = 0; i < nrows; i++)
    {
        (*A)[i] = &x[i*ncols];       
    }
}

void zeros_array(int n,double **x)
{
    int i;
    empty_array(n,x);
    for(i = 0; i < n; i++)
    {
        (*x)[i] = 0;
    }
}

void ones_array(int n,double **x)
{
    constant_array(n,x,1.0);
}

void char_array(int n, char **c)
{
    *c = malloc(n*sizeof(char));
    if (c == NULL)
    {
        print_essential("char_array : array allocation error\n");
        exit(0);
    }
}

void constant_array(int n,double **x, double value)
{
    int i;
    empty_array(n,x);
    for(i = 0; i < n; i++)
    {
        (*x)[i] = value;
    }
}

void linspace_array(double a,double b,int n,double **x)
{
    int i;
    empty_array(n,x);
    if (n == 1)
    {
        (*x)[0] = b;  /* This is what Matlab's linspace does */
        return;
    }
    double h = (b-a)/(n-1);
    for(i = 0; i < n; i++)
    {
        (*x)[i] = a + i*h;
    }
}

void random_array(int n, double **x)
{
    empty_array(n,x);    
    random_seed();   
    int i;

    for(i=0; i < n; i++)
    {
        (*x)[i] = random_number();
    }
}

void delete_array(void **x)
{
    free(*x);
}

void delete_array2(void ***A)
{
    free(**A);
    free(*A);
}



/* --------------------------------------------------------------------
    Operations on arrays

    Routines : 
        sum_array  : sum entries and return scalar.
    ----------------------------------------------------------------- */

double sum_array(int n, double *x)
{
    int i;
    double s;

    s = 0;
    for(i = 0; i < n; i++)
    {
        s += x[i];
    }
    return s;
}

/* --------------------------------------------------------------------
    Input/Output

    Read routines from the command line; 
    Print values either from node 0 only or from each processor. 

        print_global  : Print only from node 0
        print_debug   : print from each processor

    Example : 

        print_global("hello!\n");

        returns : 
        Processor [0] : hello!

        print_debug("hello!\n");

        returns : 
        Processor [0] : hello!
        Processor [1] : hello!
        Processor [2] : hello!
        Processor [3] : hello!
        
    ------------------------------------------------------------------ */

void read_int(int argc, char** argv, const char arg[], int* value,int *err)
{
    *err = 1;  /* Nothing found yet */
    int arg_index = 1;     /* Skip first argument */
    while (arg_index < argc)
    {
        if (strcmp(argv[arg_index], arg) == 0)
        {
            arg_index++;
            *value = atoi(argv[arg_index++]);   
            *err = 0;         
            return;
        }
        else
        {
            arg_index++;
        }
    }
}

void read_double(int argc, char** argv, const char arg[], double* value,int *err)
{
    *err = 1;  /* Nothing found yet */
    int arg_index = 1;     /* Skip first argument */
    while (arg_index < argc)
    {
        if (strcmp(argv[arg_index], arg) == 0)
        {
            arg_index++;
            *value = atof(argv[arg_index++]);   
            *err = 0;         
            return;
        }
        else
        {
            arg_index++;
        }
    }
}

void read_string(int argc, char** argv, const char arg[], char* value,int *err)
{
    *err = 1;  /* Nothing found yet */
    int arg_index = 1;     /* Skip first argument */
    while (arg_index < argc)
    {
        if (strcmp(argv[arg_index], arg) == 0)
        {
            arg_index++;            
            strcpy(value,argv[arg_index++]);   
            *err = 0;         
            return;
        }
        else
        {
            arg_index++;
        }
    }
}

void read_loglevel(int argc, char** argv)
{
    int err;
    char logstr[20];
    read_string(argc,argv, "--loglevel", logstr, &err);
    if (err > 0)
    {
        /* No loglevel specified at command */
        strcpy(logstr,"essential"); /* Default */
    }

    loglevel_t l;
    char loglevel_list[5][11];  /* 11 = length("production") + 1 */
    strcpy(loglevel_list[SILENT],     "silent");
    strcpy(loglevel_list[ESSENTIAL],  "essential");
    strcpy(loglevel_list[PRODUCTION], "production");
    strcpy(loglevel_list[INFO],       "info");
    strcpy(loglevel_list[DEBUG],      "debug");

    for(l = SILENT; l <= DEBUG; l++)
    {
        if (strcmp(logstr,loglevel_list[l]) == 0)    
        {
            s_loglevel = l;
            return;
        }        
    }
    s_loglevel = PRODUCTION;  /* Default, before anything else is set */
}

void print_global(const char* format, ... )
{
    if (s_rank == 0)
    {
        va_list arglist;
        printf( "Processor [0] : " );
        va_start( arglist, format );
        vprintf( format, arglist );
        va_end( arglist );
    }
}

void print_essential(const char* format, ... )
{
    /* Only print if on processor 0 */
    if (s_rank == 0 && s_loglevel > SILENT)
    {
        va_list arglist;
        printf( "Processor [0] : " );
        va_start( arglist, format );
        vprintf( format, arglist );
        va_end( arglist );
    }
}

void print_info(const char* format, ... )
{
    /* Include rank number in print statement */
    if (s_rank == 0 && s_loglevel >= INFO)
    {
        va_list arglist;
        printf( "Processor [%d] : ",s_rank);
        va_start( arglist, format );
        vprintf( format, arglist );
        va_end( arglist );
    }
}

void print_debug(const char* format, ... )
{
    /* Include rank number in print statement */
    if (s_loglevel >= DEBUG)
    {
        va_list arglist;
        printf( "Processor [%d] : ",s_rank);
        va_start( arglist, format );
        vprintf( format, arglist );
        va_end( arglist );        
    }
}

/* -------------------------------------------------
    Miscellaneous routines
   ----------------------------------------------- */ 
void set_rank(int  rank)
{
    /* This must be called so that print_debug and print_global work */
    s_rank = rank;
}

void sleep(double t_total)
{
    double t0, t1;
    t0 = clock();
    t1 = t0;
    while ((t1-t0)/CLOCKS_PER_SEC < t_total)
    {
        t1 = clock();
    }
}

double random_number()
{
  return (double) rand() / (double) RAND_MAX ;
}

/*
int random_int(int m, int n)
{
    int r = floor((n+0.999999-m)*random_number());
    return r;
}
*/

void random_seed()
{
    srand(time(NULL) + s_rank);
    int skip = random_number();  /* To start the seed?  */
}

int pow2(int p)
{
    /* Compute n = 2^p */
    int n,i;

    n = 1;
    for(i = 0; i < p; i++)
    {
        n *= 2;
    }
    return n;
}


