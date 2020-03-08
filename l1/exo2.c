#include <stdio.h>
#include <omp.h>
#include <math.h>
#include <x86intrin.h>

#define NBEXPERIMENTS    22
static long long unsigned int experiments [NBEXPERIMENTS] ;


#define N              512
#define TILE           16

typedef double vector [N] ;

typedef double matrix [N][N] ;

static vector a, b, c ;
static matrix M1, M2 ;

long long unsigned int average (long long unsigned int *exps)
{
    unsigned int i ;
    long long unsigned int s = 0 ;

    for (i = 2; i < (NBEXPERIMENTS-2); i++)
    {
        s = s + exps [i] ;
    }

    return s / (NBEXPERIMENTS-2) ;
}


void init_vector (vector X, const double val)
{
    register unsigned int i ;

    for (i = 0 ; i < N; i++)
        X [i] = val ;

    return ;
}

void init_matrix (matrix X, const double val)
{
    register unsigned int i, j;

    for (i = 0; i < N; i++)
    {
        for (j = 0 ;j < N; j++)
        {
            X [i][j] = val ;
        }
    }
}


void print_vectors (vector X, vector Y)
{
    register unsigned int i ;

    for (i = 0 ; i < N; i++)
        printf (" X [%d] = %le Y [%d] = %le\n", i, X[i], i,Y [i]) ;

    return ;
}

void add_vectors1 (vector X, vector Y, vector Z)
{
    register unsigned int i ;

#pragma omp parallel for schedule(static)
    for (i=0; i < N; i++)
        X[i] = Y[i] + Z[i];

    return ;
}

void add_vectors2 (vector X, vector Y, vector Z)
{
    register unsigned int i ;

#pragma omp parallel for schedule(dynamic)
    for (i=0; i < N; i++)
        X[i] = Y[i] + Z[i];

    return ;
}

double dot1 (vector X, vector Y)
{
    register unsigned int i ;
    register double dot ;


    dot = 0.0 ;
#pragma omp parallel for schedule(static) reduction (+:dot)
    for (i=0; i < N; i++)
        dot += X [i] * Y [i];

    return dot ;
}

double dot2 (vector X, vector Y)
{
    register unsigned int i ;
    register double dot ;


    dot = 0.0 ;
#pragma omp parallel for schedule(dynamic) reduction (+:dot)
    for (i=0; i < N; i++)
        dot += X [i] * Y [i];

    return dot ;
}

double dot3 (vector X, vector Y)
{
    register unsigned int i ;
    register double dot ;

    dot = 0.0 ;
#pragma omp parallel for schedule(static) reduction (+:dot)
    for (i = 0; i < N; i = i + 8)
    {
        dot += X [i] * Y [i];
        dot += X [i + 1] * Y [i + 1];
        dot += X [i + 2] * Y [i + 2];
        dot += X [i + 3] * Y [i + 3];

        dot += X [i + 4] * Y [i + 4];
        dot += X [i + 5] * Y [i + 5];
        dot += X [i + 6] * Y [i + 6];
        dot += X [i + 7] * Y [i + 7];
    }

    return dot ;
}

void mult_mat_vect0 (matrix M, vector b, vector c)
{
    /*
      matrix vector multiplication
      Sequential function
    */
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            c[i] += M[i][j] * b[j];
        }
    }
    return ;
}

void mult_mat_vect1 (matrix M, vector b, vector c)
{
    /*
      matrix vector multiplication
      Parallel function with static loop scheduling
    */
    for (int i = 0; i < N; ++i) {
        double total = 0.0;
#pragma omp parallel for schedule(static) reduction (+:total)
        for (int j = 0; j < N; ++j) {
            total += M[i][j] * b[j];
        }
        c[i] = total;
    }
    return ;
}

void mult_mat_vect2 (matrix M, vector b, vector c)
{
    /*
      matrix vector multiplication
      Parallel function with static loop scheduling
      unrolling internal loop
    */
    for (int i = 0; i < N; ++i) {
        double total = 0.0;
#pragma omp parallel for schedule(static) reduction (+:total)
        for (int j = 0; j < N; j += 8) {
            total += M[i][j+0] * b[j+0];
            total += M[i][j+1] * b[j+1];
            total += M[i][j+2] * b[j+2];
            total += M[i][j+3] * b[j+3];
            total += M[i][j+4] * b[j+4];
            total += M[i][j+5] * b[j+5];
            total += M[i][j+6] * b[j+6];
            total += M[i][j+7] * b[j+7];
        }
        c[i] = total;
    }
    return ;
}

void mult_mat_vect3 (matrix M, vector b, vector c)
{
    /*
      matrix vector multiplication
      Parallel function with static loop scheduling
      unrolling internal and external loops
    */
    for (int i = 0; i < N; i += 2) {
        double total0 = 0.0;
        double total1 = 0.0;
#pragma omp parallel for schedule(static) reduction (+:total0, total1)
        for (int j = 0; j < N; j += 8) {
            total0 += M[i+0][j+0] * b[j+0];
            total0 += M[i+0][j+1] * b[j+1];
            total0 += M[i+0][j+2] * b[j+2];
            total0 += M[i+0][j+3] * b[j+3];
            total0 += M[i+0][j+4] * b[j+4];
            total0 += M[i+0][j+5] * b[j+5];
            total0 += M[i+0][j+6] * b[j+6];
            total0 += M[i+0][j+7] * b[j+7];

            total1 += M[i+1][j+0] * b[j+0];
            total1 += M[i+1][j+1] * b[j+1];
            total1 += M[i+1][j+2] * b[j+2];
            total1 += M[i+1][j+3] * b[j+3];
            total1 += M[i+1][j+4] * b[j+4];
            total1 += M[i+1][j+5] * b[j+5];
            total1 += M[i+1][j+6] * b[j+6];
            total1 += M[i+1][j+7] * b[j+7];
        }
        c[i+0] = total0;
        c[i+1] = total1;
    }
    return ;
}

void mult_mat_mat0 (matrix A, matrix B, matrix C)
{
    /*
      Matrix Matrix Multiplication
      Sequential function
    */
    // A[0][0]B[0][0]+A[0][1]B[1][0]   A[0][0]B[0][1]+A[0][1]B[1][1]
    // A[1][0]B[0][0]+A[0][1]B[1][0]   A[1][0]B[0][1]+A[1][1]B[1][1]
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return ;
}


void mult_mat_mat1 (matrix A, matrix B, matrix C)
{
    /*
      Matrix Matrix Multiplication
      Parallel function with OpenMP and static scheduling
    */
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double total = 0.0;
#pragma omp parallel for schedule(static) reduction (+:total)
            for (int k = 0; k < N; k++) {
                total += A[i][k] * B[k][j];
            }
            C[i][j] = total;
        }
    }
    return ;
}

void mult_mat_mat2 (matrix A, matrix B, matrix C)
{
    /*
      Matrix Matrix Multiplication
      Parallel function with OpenMP and static scheduling
      Unrolling the inner loop
    */
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double total = 0.0;
#pragma omp parallel for schedule(static) reduction (+:total)
            for (int k = 0; k < N; k+=8) {
                total += A[i][k+0] * B[k+0][j];
                total += A[i][k+1] * B[k+1][j];
                total += A[i][k+2] * B[k+2][j];
                total += A[i][k+3] * B[k+3][j];
                total += A[i][k+4] * B[k+4][j];
                total += A[i][k+5] * B[k+5][j];
                total += A[i][k+6] * B[k+6][j];
                total += A[i][k+7] * B[k+7][j];
            }
            C[i][j] = total;
        }
    }
    return ;
}

void mult_mat_mat3 (matrix A, matrix B, matrix C)
{
    /*
      Matrix Matrix Multiplication
      Parallel function with OpenMP and static scheduling
      With tiling and unrolling
    */
    return ;
}

// currently 3209 MHz is lowest frequency ~ 3200 MHz for calculations
int main ()
{
    int nthreads, maxnthreads ;

    int tid;

    unsigned long long int start, end ;
    unsigned long long int residu ;

    unsigned long long int av ;

    double r ;

    int exp ;

    /*
       rdtsc: read the cycle counter
    */

    start = _rdtsc () ;
    end = _rdtsc () ;
    residu = end - start ;

    /*
       Sequential code executed only by the master thread
    */

    nthreads = omp_get_num_threads();
    maxnthreads = omp_get_max_threads () ;
    printf("Sequential execution: \n# threads %d\nmax threads %d \n", nthreads, maxnthreads);

    /*
      Vector Initialization
    */

    init_vector (a, 1.0) ;
    init_vector (b, 2.0) ;

    /*
      print_vectors (a, b) ;
    */


    printf ("=============== ADD ==========================================\n") ;

    init_vector (a, 1.0) ;
    init_vector (b, 2.0) ;

    /*
      print_vectors (a, b) ;
    */

    for (exp = 0 ; exp < NBEXPERIMENTS; exp++)
    {
        start = _rdtsc () ;

        add_vectors1 (c, a, b) ;

        end = _rdtsc () ;
        experiments [exp] = end - start ;
    }

    av = average (experiments) ;
    // fp/s = x, cycles/second = y, cycles/call = av-residu, fp/call = 1024, fp/s/call = ?
    // fp/call * call/cycles * cycles/second = fp/second
    // (512*2)/av-residu * y
    // 1024 * 3.2 G/av-residu = 3276.8 G / av-residu
    printf ("OpenMP static loop %Ld cycles and %f GFLOPS\n", av-residu, (3276.8)/(av-residu) ) ;

    init_vector (a, 1.0) ;
    init_vector (b, 2.0) ;

    /*
      print_vectors (a, b) ;
    */

    for (exp = 0 ; exp < NBEXPERIMENTS; exp++)
    {
        start = _rdtsc () ;

        add_vectors2 (c, a, b) ;

        end = _rdtsc () ;
        experiments [exp] = end - start ;
    }

    av = average (experiments) ;

    printf ("OpenMP dynamic loop %Ld cycles and %f GFLOPS\n", av-residu, (3276.8)/(av-residu) ) ;

    printf ("==============================================================\n") ;

    printf ("====================DOT =====================================\n") ;

    init_vector (a, 1.0) ;
    init_vector (b, 2.0) ;

    for (exp = 0 ; exp < NBEXPERIMENTS; exp++)
    {
        start = _rdtsc () ;

        r = dot1 (a, b) ;

        end = _rdtsc () ;
        experiments [exp] = end - start ;
    }

    av = average (experiments) ;
    // fp/call * call/cycles * cycles/second = fp/second
    // (512*3)/av-residu * y
    // 1536 * 3.2 G/av-residu = 4915.2 G / av-residu
    printf ("OpenMP static loop dot %e: %Ld cycles and %f GFLOPS\n", r, av-residu, (4915.2)/(av-residu)) ;

    init_vector (a, 1.0) ;
    init_vector (b, 2.0) ;

    for (exp = 0 ; exp < NBEXPERIMENTS; exp++)
    {
        start = _rdtsc () ;

        r = dot2 (a, b) ;

        end = _rdtsc () ;
        experiments [exp] = end - start ;
    }

    av = average (experiments) ;

    printf ("OpenMP dynamic loop dot %e: %Ld cycles and %f GFLOPS\n", r, av-residu, (4915.2)/(av-residu)) ;

    init_vector (a, 1.0) ;
    init_vector (b, 2.0) ;

    for (exp = 0 ; exp < NBEXPERIMENTS; exp++)
    {
        start = _rdtsc () ;

        r = dot3 (a, b) ;

        end = _rdtsc () ;
        experiments [exp] = end - start ;
    }

    av = average (experiments) ;
    //
    printf ("OpenMP static unrolled loop dot %e: %Ld cycles and %f GFLOPS\n", r, av-residu, (4915.2)/(av-residu)) ;

    printf ("=============================================================\n") ;

    printf ("======================== Mult Mat Vector =====================================\n") ;

    init_matrix (M1, 1.0) ;
    init_vector (b, 2.0) ;

    for (exp = 0 ; exp < NBEXPERIMENTS; exp++)
    {
        start = _rdtsc () ;

        mult_mat_vect0 (M1, b, a) ;

        end = _rdtsc () ;
        experiments [exp] = end - start ;
    }

    av = average (experiments) ;
    // fp/call * call/cycles * cycles/second = fp/second
    // (512*512*3)/av-residu * y
    // 786432 * 3.2 G/av-residu = 2516582.4 M / av-residu
    printf ("OpenMP static loop MultMatVect0: %Ld cycles and %f MFLOPS\n", av-residu, (2516582.4)/(av-residu)) ;

    init_matrix (M1, 1.0) ;
    init_vector (b, 2.0) ;

    for (exp = 0 ; exp < NBEXPERIMENTS; exp++)
    {
        start = _rdtsc () ;

        mult_mat_vect1 (M1, b, a) ;

        end = _rdtsc () ;
        experiments [exp] = end - start ;
    }

    av = average (experiments) ;

    printf ("OpenMP static loop MultMatVect1: %Ld cycles and %f MFLOPS\n", av-residu, (2516582.4)/(av-residu)) ;

    init_matrix (M1, 1.0) ;
    init_vector (b, 2.0) ;

    for (exp = 0 ; exp < NBEXPERIMENTS; exp++)
    {
        start = _rdtsc () ;

        mult_mat_vect2 (M1, b, a) ;

        end = _rdtsc () ;
        experiments [exp] = end - start ;
    }

    av = average (experiments) ;

    printf ("OpenMP static loop MultMatVect2: %Ld cycles and %f MFLOPS\n", av-residu, (2516582.4)/(av-residu)) ;

    init_matrix (M1, 1.0) ;
    init_vector (b, 2.0) ;

    for (exp = 0 ; exp < NBEXPERIMENTS; exp++)
    {
        start = _rdtsc () ;

        mult_mat_vect3 (M1, b, a) ;

        end = _rdtsc () ;
        experiments [exp] = end - start ;
    }

    av = average (experiments) ;

    printf ("OpenMP static loop MultMatVect3: %Ld cycles and %f MFLOPS\n", av-residu, (2516582.4)/(av-residu)) ;

    printf ("==============================================================\n") ;

    printf ("======================== Mult Mat Mat =====================================\n") ;

    init_matrix (M1, 1.0) ;
    init_matrix (M2, 2.0) ;

    for (exp = 0 ; exp < NBEXPERIMENTS; exp++)
    {
        start = _rdtsc () ;

        mult_mat_mat0 (M1, M2, M2) ;

        end = _rdtsc () ;
        experiments [exp] = end - start ;
    }

    av = average (experiments) ;
    // fp/call * call/cycles * cycles/second = fp/second
    // (512*512*512*3)/av-residu * y
    // 402623184 * 3.2 G/av-residu = 1288490188.8 M / av-residu
    printf ("Sequential matrice vector multiplication:\t %Ld cycles and %f MFLOPS\n", av-residu, (1288490188.8)/(av-residu)) ;

    init_matrix (M1, 1.0) ;
    init_matrix (M2, 2.0) ;


    for (exp = 0 ; exp < NBEXPERIMENTS; exp++)
    {
        start = _rdtsc () ;

        mult_mat_mat1 (M1, M2, M2) ;

        end = _rdtsc () ;
        experiments [exp] = end - start ;
    }

    av = average (experiments) ;

    printf ("OpenMP static loop MultMatVect1:\t\t %Ld cycles and %f MFLOPS\n", av-residu, (1288490188.8)/(av-residu)) ;

    init_matrix (M1, 1.0) ;
    init_matrix (M2, 2.0) ;

    for (exp = 0 ; exp < NBEXPERIMENTS; exp++)
    {
        start = _rdtsc () ;

        mult_mat_mat2 (M1, M2, M2) ;

        end = _rdtsc () ;
        experiments [exp] = end - start ;
    }

    av = average (experiments) ;


    printf ("OpenMP unrolled loop MultMatMat2:\t\t %Ld cycles and %f MFLOPS\n", av-residu, (1288490188.8)/(av-residu)) ;

    init_matrix (M1, 1.0) ;
    init_matrix (M2, 2.0) ;

    for (exp = 0 ; exp < NBEXPERIMENTS; exp++)
    {
        start = _rdtsc () ;

        mult_mat_mat3 (M1, M2, M2) ;

        end = _rdtsc () ;
        experiments [exp] = end - start ;
    }

    av = average (experiments) ;

    printf ("OpenMP Tiled loop MultMatMat3:\t\t\t %Ld cycles and %f MFLOPS\n", av-residu, (1288490188.8)/(av-residu)) ;

    return 0;

}

