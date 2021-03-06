#include <stdio.h>
#include <omp.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

#include <x86intrin.h>

#include "sorting.h"

static int compare (const void *x, const void *y)
{
	/* cast x and y to uint64_t* before comparing */
	uint64_t xVal = *((uint64_t*)x);
	uint64_t yVal = *((uint64_t*)y);
	return xVal - yVal;
}

void do_quick_sort(uint64_t *T, const int size){
	qsort(T, size, sizeof(uint64_t), compare);
}

void sequential_qsort_sort (uint64_t *T, const int size)
{
	do_quick_sort(T, size);
	return ;
}

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



void parallel_qsort_sort (uint64_t *T, const uint64_t size, const uint64_t blkSiz)
{
	int nbUnit = size/blkSiz;
	uint64_t blockSiz = blkSiz;
#pragma omp parallel for schedule(static)
	for(int i = 0; i < nbUnit; ++i){
		uint64_t *cur = &T[i*blockSiz];
		do_quick_sort(cur, blockSiz);
	}
	while(nbUnit != 0){
		nbUnit = nbUnit / 2;
		for(int i = 0; i < nbUnit; ++i){
			uint64_t *cur = &T[i*(blockSiz*2)];
			merge(cur, blockSiz);
		}
		blockSiz = blockSiz * 2;
	}
	return;
}


void parallel_qsort_sort1 (uint64_t *T, const uint64_t size, const uint64_t blkSiz)
{
	int nbUnit = size/blkSiz;
	uint64_t blockSiz = blkSiz;
#pragma omp parallel for schedule(static)
	for(int i = 0; i < nbUnit; ++i){
		uint64_t *cur = &T[i*blockSiz];
		do_quick_sort(cur, blockSiz);
	}
	while(nbUnit != 0){
		nbUnit = nbUnit / 2;
#pragma omp parallel for schedule(static)
		for(int i = 0; i < nbUnit; ++i){
			uint64_t *cur = &T[i*(blockSiz*2)];
			merge(cur, blockSiz);
		}
		blockSiz = blockSiz * 2;
	}
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
		fprintf (stderr, "qsort.run N(lg(array size)) L(lg(block size)) \n") ;
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

	for (exp = 0 ; exp < NBEXPERIMENTS; exp++){
#ifdef RINIT
		init_array_random (X, sN);
#else
		init_array_sequence (X, sN);
#endif


		start = _rdtsc () ;

		sequential_qsort_sort (X, sN) ;

		end = _rdtsc () ;
		experiments [exp] = end - start ;

		/* verifying that X is properly sorted */
#ifdef RINIT
		if (! is_sorted (X, sN))
#else
		if (! is_sorted_sequence (X, sN))
#endif
		{
			fprintf(stderr, "ERROR: the sequential sorting of the array failed\n") ;
			print_array (X, sN) ;
			exit (-1) ;
		}
	}

	av = average_time() ;

	printf ("\n qsort serial \t\t\t %.2lf Mcycles\n\n", (double)av/1000000) ;


	for (exp = 0 ; exp < NBEXPERIMENTS; exp++)
	{
#ifdef RINIT
		init_array_random (X, N);
#else
		init_array_sequence (X, N);
#endif

		start = _rdtsc () ;

		parallel_qsort_sort (X, N, L) ;


		end = _rdtsc () ;
		experiments [exp] = end - start ;

		/* verifying that X is properly sorted */
#ifdef RINIT
		if (! is_sorted (X, N))
#else
		if (! is_sorted_sequence (X, N))
#endif
		{
			fprintf(stderr, "ERROR: the parallel sorting(with sequential merge) of the array failed\n") ;
			print_array (X, sN) ;
			exit (-1) ;
		}


	}

	av = average_time() ;
	printf ("\n qsort parallel (merge seq) \t %.2lf Mcycles\n\n", (double)av/1000000) ;

	for (exp = 0 ; exp < NBEXPERIMENTS; exp++)
	{
#ifdef RINIT
		init_array_random (X, N);
#else
		init_array_sequence (X, N);
#endif

		start = _rdtsc () ;

		parallel_qsort_sort1 (X, N, L) ;

		end = _rdtsc () ;
		experiments [exp] = end - start ;

		/* verifying that X is properly sorted */
#ifdef RINIT
		if (! is_sorted (X, N))
#else
		if (! is_sorted_sequence (X, N))
#endif
		{
			fprintf(stderr, "ERROR: the parallel sorting(with parallel merge) of the array failed\n") ;
			exit (-1) ;
		}


	}

	av = average_time() ;
	printf ("\n qsort parallel \t\t %.2lf Mcycles\n\n", (double)av/1000000) ;

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

	sequential_qsort_sort (Y, N) ;
	parallel_qsort_sort1 (Z, N, L) ;

	if (! are_vector_equals (Y, Z, N)) {
		fprintf(stderr, "ERROR: sorting with the sequential and the parallel algorithm does not give the same result\n") ;
		exit (-1) ;
	}


	free(X);
	free(Y);
	free(Z);
#endif
}
