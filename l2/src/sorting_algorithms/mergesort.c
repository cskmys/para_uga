#include <stdio.h>
#include <omp.h>
#include <stdint.h>
#include <string.h>

#include <x86intrin.h>

#include "sorting.h"

/* 
   Merge two sorted chunks of array T!
   The two chunks are of size size
   First chunck starts at T[0], second chunck starts at T[size]
 */
void merge (uint64_t *T, const uint64_t size)
{
	uint64_t *X = (uint64_t *) malloc (2 * size * sizeof(uint64_t)) ;

	uint64_t i = 0 ;
	uint64_t j = size ;
	uint64_t k = 0 ;

	while ((i < size) && (j < 2*size))
	{
		if (T[i] < T [j])
		{
			X [k] = T [i] ;
			i = i + 1 ;
		}
		else
		{
			X [k] = T [j] ;
			j = j + 1 ;
		}
		k = k + 1 ;
	}

	if (i < size)
	{
		for (; i < size; i++, k++)
		{
			X [k] = T [i] ;
		}
	}
	else
	{
		for (; j < 2*size; j++, k++)
		{
			X [k] = T [j] ;
		}
	}

	memcpy (T, X, 2*size*sizeof(uint64_t)) ;
	free (X) ;

	return ;
}
/*
   merge sort -- sequential, parallel -- 
 */

void sequential_merge_sort (uint64_t *T, const uint64_t size)
{
	/* TODO: sequential implementation of merge sort */

	return ;
}

void parallel_merge_sort (uint64_t *T, const uint64_t size, const uint64_t blkSiz)
{
	/* TODO: parallel implementation of merge sort */

	return;
}


int main (int argc, char **argv)
{
	uint64_t start, end;
	uint64_t av ;
	unsigned int exp ;

	/* the program takes one parameter N which is the size of the array to
       be sorted. The array will have size 2^N */
	if (argc != 3)
	{
		fprintf (stderr, "mergesort.run N(lg(array size)) L(lg(block size)) \n") ;
		exit (-1) ;
	}

	uint64_t N = 1 << (atoi(argv[1])) ;
	uint64_t L = 1 << (atoi(argv[2])) ; // to keep code simple we assume both params are powers of 2

	if(N <= L){
		fprintf (stderr, "provide an array size > block size \n") ;
		exit (-1) ;
	}
	/* the array to be sorted */
	uint64_t *X = (uint64_t *) malloc (N * sizeof(uint64_t)) ;

	printf("--> Sorting an array of size %lu\n",N);

	uint64_t sN = N;

#ifdef TEST_ONLY_PARA
	sN = TEST_ONLY_PARA_N;
	printf("--> Sorting an array of size %lu sequentially and array of size %lu in parallel\n", sN, N);
#else
	printf("--> Sorting an array of size %lu\n",N);
#endif
#ifdef RINIT
	printf("--> The array is initialized randomly\n");
#endif
	for (exp = 0 ; exp < NBEXPERIMENTS; exp++)
	{
#ifdef RINIT
		init_array_random (X, sN);
#else
		init_array_sequence (X, sN);
#endif


		start = _rdtsc () ;

		sequential_merge_sort (X, sN) ;

		end = _rdtsc () ;
		experiments [exp] = end - start ;

		/* verifying that X is properly sorted */
#ifdef RINIT
		if (! is_sorted (X, sN))
		{
			fprintf(stderr, "ERROR: the sequential sorting of the array failed\n") ;
			print_array (X, sN) ;
			exit (-1) ;
		}
#else
		if (! is_sorted_sequence (X, sN))
		{
			fprintf(stderr, "ERROR: the sequential sorting of the array failed\n") ;
			print_array (X, sN) ;
			exit (-1) ;
		}
#endif
	}

	av = average_time() ;

	printf ("\n mergesort serial \t\t\t %.2lf Mcycles\n\n", (double)av/1000000) ;


	for (exp = 0 ; exp < NBEXPERIMENTS; exp++)
	{
#ifdef RINIT
		init_array_random (X, N);
#else
		init_array_sequence (X, N);
#endif

		start = _rdtsc () ;

		parallel_merge_sort (X, N, L) ;

		end = _rdtsc () ;
		experiments [exp] = end - start ;

		/* verifying that X is properly sorted */
#ifdef RINIT
		if (! is_sorted (X, N))
		{
			fprintf(stderr, "ERROR: the parallel sorting of the array failed\n") ;
			exit (-1) ;
		}
#else
		if (! is_sorted_sequence (X, N))
		{
			fprintf(stderr, "ERROR: the parallel sorting of the array failed\n") ;
			exit (-1) ;
		}
#endif


	}

	av = average_time() ;
	printf ("\n mergesort parallel \t\t %.2lf Mcycles\n\n", (double)av/1000000) ;

	/* print_array (X, N) ; */
#ifndef TEST_ONLY_PARA
	/* before terminating, we run one extra test of the algorithm */
	uint64_t *Y = (uint64_t *) malloc (N * sizeof(uint64_t)) ;
	uint64_t *Z = (uint64_t *) malloc (N * sizeof(uint64_t)) ;

#ifdef RINIT
	init_array_random (Y, N);
#else
	init_array_sequence (Y, N);
#endif

	memcpy(Z, Y, N * sizeof(uint64_t));

	sequential_merge_sort (Y, N) ;
	parallel_merge_sort (Z, N, L) ;

	if (! are_vector_equals (Y, Z, N)) {
		fprintf(stderr, "ERROR: sorting with the sequential and the parallel algorithm does not give the same result\n") ;
		print_array (Y, N) ;
		print_array (Z, N) ;
		exit (-1) ;
	}


	free(X);
	free(Y);
	free(Z);
#endif
}
