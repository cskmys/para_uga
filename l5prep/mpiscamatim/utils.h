#ifndef UTILS_H
#define UTILS_H

#include <math.h>
#include <assert.h>

#define MAX_NB_PROC 64

#define NBEXPERIMENTS    5

typedef double mat;
#define PRINT_FMT_SPECIFIER "%.2f "

#define NULL_PTR_CHK(x) assert(x != NULL)
#define SAME_PTR_CHK(x, y) assert(x != y)
#define ABOVE_ZERO_CHK(x) assert(x > 0.0)

#define ELE_SIZ(x) sizeof(x[0])
#define NB_ELE(x) sizeof(x)/ELE_SIZ(x)

void init_comm_pipe(int rank);
void deinit_comm_pipe(int rank);

void mpi_printf(const char *fmt, ...);
void mpi_scr_printf(const char *fmt, ...);

/* return the average time in cycles over the values stored in
 * experiments vector */
double average_time( void );
void set_exp_time(int trial, double tim);

#endif // UTILS_H