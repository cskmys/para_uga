#ifndef FOX_H
#define FOX_H

#include "utils.h"

void fox_matrix_mul(int size, int rank, mat *A, int rA, int cA, mat *B, int rB, int cB, mat *C);

void fox_mat_mul_prep(int size, int rank, mat *A, int rA, int cA, mat *B, int rB, int cB, mat *C);

#endif