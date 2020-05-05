#ifndef MATRIX_H
#define MATRIX_H

#include "utils.h"
#include <stdlib.h>

void chkMatrixValidity(mat *matr, size_t row, size_t col);
void printMatrix(mat *buf, size_t row, size_t col);
mat *allocMatMem(size_t r, size_t c);
void memsetMatrix(mat *matr, int val, size_t row, size_t col);
void randsetMatrix(mat *matr, int maxVal, size_t row, size_t col);
void setIdentityMatrix(mat *matr, size_t dim);
void cpyMatrix(mat *dst, mat *src, size_t row, size_t col, size_t cpyStartRowIdx, size_t nbRowsCpy);
void sumMatrix(mat *x, mat *h, size_t r, size_t c, mat *y);
// void matrixMul(mat *matA, size_t rowA, size_t colA, mat *matB, size_t rowB, size_t colB, mat **matC, size_t *rowC, size_t *colC);
void matrixMul(mat *matA, size_t rowA, size_t colA, mat *matB, size_t rowB, size_t colB, mat *matC);
#endif // MATRIX_H
