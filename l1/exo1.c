#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <sched.h>
#include <omp.h>

#define NB_THREADS 16

void Hello(char *st)
{

	int my_rank, nthreads, nmaxthreads;

	/*
      to get the max number of threads
	 */

	nmaxthreads = omp_get_max_threads () ;

	/*
     to get the current number of threads
	 */
	nthreads = omp_get_num_threads();

	/*
      to get my rank
	 */
	my_rank = omp_get_thread_num();

	printf("%s thread %d running on cpu %d of team %d (max_num_threads is %d)\n", st, my_rank, sched_getcpu(), nthreads, nmaxthreads);

}


int main (int argc, char*argv[])
{

	printf("I can execute a max of %d threads in parallel\n", omp_get_num_procs());
	/*
      Program starts here with the master thread
	 */

	int nthreads = NB_THREADS;
	printf("I am the master thread %d running on cpu %d and I start\n", sched_getcpu(), omp_get_thread_num ());

	/*
       This is a parallel region
       All threads will execute this region
	 */
	printf ("Starting Region 1 \n") ;
#pragma omp parallel num_threads (nthreads)
	{
		Hello("Region 1 ") ;
	}

	printf ("End of Region 1\n") ;


	printf ("Starting Region 2 \n") ;


	/*
     This is a parallel region
     half of threads will execute this region
	 */
#pragma omp parallel num_threads (nthreads/2)
	{
		Hello ("Region 2 ") ;
	}
	/*
      End of the parallel region
      The master thread is alone in this serial region
	 */
	printf ("End of Region 2\n") ;


	/*
       This is a parallel region
       quarter of threads will execute this region
	 */
	printf ("Starting Region 3 \n") ;


#pragma omp parallel num_threads (nthreads/4)
	{
		Hello ("Region 3 ") ;
	}

	/*
      End of the parallel region
      The master thread is alone in this serial region
	 */
	printf ("End of Region 3\n") ;

	printf("I am the master thread %d running on cpu %d and I complete\n", sched_getcpu(), omp_get_thread_num ());

	return 0;
}
