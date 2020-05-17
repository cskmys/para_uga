#ifndef RING_H
#define RING_H

#include "utils.h"

void ring_matrix_mul(int size, int rank, mat *A, int rA, int cA, mat *B, int rB, int cB, mat *C);

#endif