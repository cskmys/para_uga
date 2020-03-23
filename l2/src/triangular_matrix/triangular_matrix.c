#include <stdio.h>
#include <math.h>
#include <omp.h>

#include <x86intrin.h>

#define NBEXPERIMENTS    7

static long long unsigned int experiments [NBEXPERIMENTS] ;

#define N              1024
#define CHUNK             16
#define NB_OF_THREADS     8


typedef double vector [N] ;

typedef double matrix [N][N] ;

matrix M ;
vector v1 ;
vector v2 ;

long long unsigned int average (long long unsigned int *exps)
{
	unsigned int i ;
	long long unsigned int s = 0 ;

	for (i = 2; i < (NBEXPERIMENTS-2); i++)
	{
		s = s + exps [i] ;
	}

	return s / (NBEXPERIMENTS-4) ;
}

void init_vector (vector X, const double val)
{
	register unsigned int i ;

	for (i = 0 ; i < N; i++)
		X [i] = val ;

	return ;
}

void init_matrix_inf (matrix X, double val)
{
	register unsigned int i, j;

	for (i = 0; i < N ; i++)
	{
		for (j = 0 ; j < i; j++)
		{
			X [i][j] = val ;
			X [j][i] = 0.0 ;
		}
		X [i][i] = val ;
	}
}

void print_vector (vector X)
{
	register unsigned int i ;

	for (i = 0 ; i < N; i++)
		printf (" %3.2f", X[i]) ;

	printf ("\n\n") ;
	fflush (stdout) ;

	return ;
}

void print_matrix (matrix Mat)
{
	register unsigned int i, j ;

	for (i = 0 ; i < N; i++)
	{
		for (j = 0 ; j < N; j++)
		{
			printf (" %3.2f ", Mat[i][j]) ;
		}
		printf ("\n") ;
	}
	printf ("\n") ;
	return ;
}



void mult_mat_vector (matrix Mat, vector b, vector c)
{
	register unsigned int i ;
	register unsigned int j ;
	register double r ;

	for ( i = 0 ; i < N ; i = i + 1)
	{
		r = 0.0 ;
		for (j = 0 ; j < N ; j = j + 1)
		{
			r += Mat [i][j] * b [j] ;
		}
		c [i] = r ;
	}
	return ;
}

void mult_mat_vector_tri_inf (matrix Mat, vector b, vector c)
{
	/*
      this function is sequential (no OpenMP directive)
      Computes the Multiplication between the vector b and the Triangular Lower Matrix
	 */
	for (register int i = 0; i < N; ++i) {
		register double r = 0.0;
		for (register int j = 0; j <= i; ++j) {
			r += Mat[i][j] * b[j];
		}
		c[i] = r;
	}
	return ;
}

void mult_mat_vector_tri_inf1 (matrix Mat, vector b, vector c)
{
	/*
      this function is parallel (with OpenMP directive, static scheduling)
       Computes the Multiplication between the vector b and the Triangular Lower Matrix
	 */
	register double total=0.0;
#pragma omp parallel for schedule(static) reduction (+:total) num_threads(NB_OF_THREADS)
	for (register int i = 0; i < N; ++i) {
		 total = 0.0;
		for (register int j = 0; j <= i; ++j) {
			total += Mat[i][j] * b[j];
		}
		c[i] = total;
	}
	return ;
}

void mult_mat_vector_tri_inf2 (matrix Mat, vector b, vector c)
{
	/*
      this function is parallel (with OpenMP directive, dynamic scheduling)
       Computes the Multiplication between the vector b and the Triangular Lower Matrix
	 */
	register double total=0.0;
#pragma omp parallel for schedule(dynamic,CHUNK) reduction (+:total) num_threads(2)
	for (register int i = 0; i < N; ++i) {
		 total = 0.0;
		for (register int j = 0; j <= i; ++j) {
			total += Mat[i][j] * b[j];
		}
		c[i] = total;
	}
	return ;
}

void mult_mat_vector_tri_inf3 (matrix Mat, vector b, vector c)
{
	/*
      this function is parallel (with OpenMP directive, guided scheduling)
      Computes the Multiplication between the vector b and the Triangular Lower Matrix
	 */
	register double total=0.0;
#pragma omp parallel for schedule(guided) reduction (+:total)  num_threads(2)
	for (register int i = 0; i < N; ++i) {
		 total = 0.0;
		for (register int j = 0; j <= i; ++j) {
			total += Mat[i][j] * b[j];
		}
		c[i] = total;
	}
	return ;
}

void mult_mat_vector_tri_inf4 (matrix Mat, vector b, vector c)
{
	/*
      this function is parallel (with OpenMP directive, runtime scheduling)
      Computes the Multiplication between the vector b and the Triangular Lower Matrix
	 */
	register double total=0.0;
#pragma omp parallel for schedule(runtime) reduction (+:total) num_threads(NB_OF_THREADS)
	for (register int i = 0; i < N; ++i) {
		 total = 0.0;
		for (register int j = N - 1 - i; j < N; ++j) {
			total += Mat[i][j] * b[j];
		}
		c[i] = total;
	}
	return ;
}

int main ( void )
{
	unsigned long long start, end ;
	unsigned long long residu ;

	unsigned long long av ;
	unsigned int exp ;

	printf ("number of threads %d\n", NB_OF_THREADS) ;

	/* rdtsc: read the cycle counter */
	start = _rdtsc () ;
	end = _rdtsc () ;
	residu = end - start ;

	init_vector (v1, 1.0) ;
	init_matrix_inf (M, 2.0) ;

	/*
      print_vector (v1) ;
      print_matrix (M) ;
	 */

	for (exp = 0 ; exp < NBEXPERIMENTS; exp++)
	{
		start = _rdtsc () ;

		mult_mat_vector (M, v1, v2)   ;

		end = _rdtsc () ;
		experiments [exp] = end - start ;
	}

	av = average (experiments) ;

	/*
      print_vector (v2) ;
	 */

	printf ("Full matrix multiplication vector \t\t  %Ld cycles\n", av-residu) ;

	init_vector (v1, 1.0) ;
	init_matrix_inf (M, 2.0) ;


	for (exp = 0 ; exp < NBEXPERIMENTS; exp++)
	{
		start = _rdtsc () ;

		mult_mat_vector_tri_inf (M, v1, v2)  ;

		end = _rdtsc () ;
		experiments [exp] = end - start ;
	}

	av = average (experiments) ;

	/*
      print_vector (v2) ;
	 */

	printf ("Triangular Matrix multiplication vector\t\t  %Ld cycles\n", av-residu) ;

	init_vector (v1, 1.0) ;
	init_matrix_inf (M, 2.0) ;

	for (exp = 0 ; exp < NBEXPERIMENTS; exp++)
	{
		start = _rdtsc () ;

		mult_mat_vector_tri_inf1 (M, v1, v2)  ;

		end = _rdtsc () ;
		experiments [exp] = end - start ;
	}

	av = average (experiments) ;

	/*
      print_vector (a) ;
	 */

	printf ("Parallel Loop Static Scheduling \t\t  %Ld cycles\n", av-residu) ;

	init_vector (v1, 1.0) ;
	init_matrix_inf (M, 2.0) ;

	for (exp = 0 ; exp < NBEXPERIMENTS; exp++)
	{
		start = _rdtsc () ;

		mult_mat_vector_tri_inf2 (M, v1, v2)  ;

		end = _rdtsc () ;
		experiments [exp] = end - start ;
	}

	av = average (experiments) ;

	/*
      print_vector (a) ;
	 */

	printf ("Parallel Loop Dynamic Scheduling \t\t  %Ld cycles\n", av-residu) ;

	init_vector (v1, 1.0) ;
	init_matrix_inf (M, 2.0) ;

	for (exp = 0 ; exp < NBEXPERIMENTS; exp++)
	{
		start = _rdtsc () ;

		mult_mat_vector_tri_inf3 (M, v1, v2)  ;

		end = _rdtsc () ;
		experiments [exp] = end - start ;
	}

	av = average (experiments) ;

	/*
       print_vector (v2) ;
	 */

	printf ("Parallel Loop Guided Scheduling \t\t  %Ld cycles\n", av-residu) ;

	init_vector (v1, 1.0) ;
	init_matrix_inf (M, 2.0) ;

	for (exp = 0 ; exp < NBEXPERIMENTS; exp++)
	{
		start = _rdtsc () ;

		mult_mat_vector_tri_inf4 (M, v1, v2)  ;

		end = _rdtsc () ;
		experiments [exp] = end - start ;
	}

	av = average (experiments) ;

	/*
       print_vector (v2) ;
	 */

	printf ("Parallel Loop Runtime Scheduling \t\t  %Ld cycles\n", av-residu) ;

	return 0;

}

