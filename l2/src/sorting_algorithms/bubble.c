#include <stdio.h>
#include <omp.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>
#include <x86intrin.h>

#include "sorting.h"

bool comp_and_swap(uint64_t *T1, uint64_t *T2){
	bool ret = false;
	if(*T1 > *T2){
		uint64_t temp = *T1;
		*T1 = *T2;
		*T2 = temp;
		ret = true;
	}
#ifdef PRINT_INT
	printf("%lu %lu %d\n", *T1, *T2, ret);
#endif
	return ret;
}

bool do_one_bubble( uint64_t *T, const uint64_t size ){
#ifdef PRINT_INT
	print_array (T, size) ;
#endif
	bool sorted = true;
	for (int i = 0; i < size - 1; ++i) {
		sorted &= !comp_and_swap(&T[i], &T[i+1]);
	}
#ifdef PRINT_INT
	print_array (T, size) ;
#endif
	return sorted;
}
/* 
   bubble sort -- sequential, parallel -- 
 */
void sequential_bubble_sort (uint64_t *T, const uint64_t size)
{
	while(do_one_bubble(T, size) == false);
	return ;
}

void parallel_bubble_sort (uint64_t *T, const uint64_t size, const uint64_t blkSiz)
{
	int sorted = 0;
	int nbUnit = size/blkSiz;
	while(sorted == 0){
		sorted = 1;
#pragma omp parallel for schedule(static) reduction (*:sorted)
		for(int i = 0; i < nbUnit; ++i){
			uint64_t *cur = &T[i*blkSiz];
			sorted *= (int) do_one_bubble(cur, blkSiz);
		}
#pragma omp parallel for schedule(static) reduction (*:sorted)
		for(int i = 0; i < nbUnit - 1; ++i){
			uint64_t *nxtBlkStart = &T[(i+1)*blkSiz];
			uint64_t *prevBlkEnd = nxtBlkStart - 1;
			sorted *= (int)!comp_and_swap(prevBlkEnd, nxtBlkStart);
		}
	}
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
		fprintf (stderr, "bubble.run N(lg(array size)) L(lg(block size)) \n") ;
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
	for (exp = 0 ; exp < NBEXPERIMENTS; exp++)
	{
#ifdef RINIT
		init_array_random (X, sN);
#else
		init_array_sequence (X, sN);
#endif


		start = _rdtsc () ;

		sequential_bubble_sort (X, sN) ;

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

	printf ("\n bubble serial \t\t\t %.2lf Mcycles\n\n", (double)av/1000000) ;

	for (exp = 0 ; exp < NBEXPERIMENTS; exp++)
	{
#ifdef RINIT
		init_array_random (X, N);
#else
		init_array_sequence (X, N);
#endif

		start = _rdtsc () ;

		parallel_bubble_sort (X, N, L) ;

		end = _rdtsc () ;
		experiments [exp] = end - start ;

		/* verifying that X is properly sorted */
#ifdef RINIT
		if (! is_sorted (X, N))
#else
		if (! is_sorted_sequence (X, N))
#endif
		{
			fprintf(stderr, "ERROR: the parallel sorting of the array failed\n") ;
			print_array (X, N) ;
			exit (-1) ;
		}

	}

	av = average_time() ;
	printf ("\n bubble parallel \t\t %.2lf Mcycles\n\n", (double)av/1000000) ;

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

	sequential_bubble_sort (Y, N) ;
	parallel_bubble_sort (Z, N, L) ;

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
