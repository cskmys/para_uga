#include "matrix_op.h"
#include <assert.h>
#include <stdio.h>
#include <string.h>

void chkMatrixValidity(mat *matr, size_t row, size_t col){
    NULL_PTR_CHK(matr);
    ABOVE_ZERO_CHK(row);
    ABOVE_ZERO_CHK(col);
}

void printMatrix(mat *buf, size_t row, size_t col){
    chkMatrixValidity(buf, row, col);

    for (size_t i = 0; i < row; ++i) {
        for (size_t j = 0; j < col; ++j) {
            printf("%d ", buf[(i * col) + j]);
        }
        printf("\n");
    }
    printf("\n");
}

mat *allocMatMem(size_t r, size_t c){
    ABOVE_ZERO_CHK(r);
    ABOVE_ZERO_CHK(c);

    return malloc(r * c * sizeof (mat));
}

void memsetMatrix(mat *matr, int val, size_t row, size_t col){
    chkMatrixValidity(matr, row, col);
    for (size_t i = 0; i < row; ++i) {
        for (size_t j = 0; j < col; ++j) {
            matr[(i * col) + j] = val;
        }
    }
}

void setIdentityMatrix(mat *matr, size_t dim){
    chkMatrixValidity(matr, dim, dim);
    for (size_t i = 0; i < dim; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            matr[(i * dim) + j] = 0;
            if(i == j){
                matr[(i * dim) + j] = 1;
            }
        }
    }
}

void cpyMatrix(mat *dst, mat *src, size_t srcRow, size_t srcCol, size_t cpyStartRowIdx, size_t nbRowsCpy){
    chkMatrixValidity(src, srcRow, srcCol);
    NULL_PTR_CHK(dst);
    SAME_PTR_CHK(dst, src);
    ABOVE_ZERO_CHK(nbRowsCpy);
    assert(!(cpyStartRowIdx < 0));
    assert(cpyStartRowIdx + nbRowsCpy <= srcRow);
    int offset = cpyStartRowIdx * srcCol;
    for (size_t i = 0; i < nbRowsCpy; ++i) {
        for (size_t j = 0; j < srcCol; ++j) {
            dst[(i * srcCol) + j] = src[(i * srcCol) + j + offset];
        }
    }
}

void sumMatrix(mat *x, mat *h, size_t r, size_t c, mat *y){
    chkMatrixValidity(x, r, c);
    NULL_PTR_CHK(h);
    NULL_PTR_CHK(y);

    for (size_t i = 0; i < r; ++i) {
        for (size_t j = 0; j < c; ++j) {
            y[((i * c) + j)] = x[(i * c) + j] + h[(i * c) + j];
        }
    }
}

void getMatAEle(mat *matr, size_t col, size_t rowNo, mat *arr){
    for (size_t j = 0; j < col; ++j) {
        arr[j] = matr[ (rowNo * col) + j];
    }
}

void getMatBEle(mat *matr, size_t nbRow, size_t nbCol, size_t colNo, mat *arr){
    for (size_t i = 0; i < nbRow; ++i) {
        arr[i] = matr[ (i * nbCol) + colNo];
    }
}

int multEle(mat *matA, mat *matB, size_t nbEle){
    int sum = 0;
    for (size_t i = 0; i < nbEle; ++i) {
        sum += (matA[i] * matB[i]);
    }
    return sum;
}

void matrixMul(mat *matA, size_t rowA, size_t colA, mat *matB, size_t rowB, size_t colB, mat **matC, size_t *rowC, size_t *colC){
    chkMatrixValidity(matA, rowA, colA);
    chkMatrixValidity(matB, rowB, colB);
    NULL_PTR_CHK(matC);
    NULL_PTR_CHK(rowC);
    NULL_PTR_CHK(colC);
    SAME_PTR_CHK(*matC,(void*)matA);
    SAME_PTR_CHK(*matC,(void*)matB);
    assert(colA == rowB);

    *rowC = rowA;
    *colC = colB;
    *matC = (void *)allocMatMem(rowA, colB);

    mat *p = *matC;

    for (size_t i = 0; i < rowA; ++i) {
        for (size_t j = 0; j < colB; ++j) {
            mat *arrA = allocMatMem(1, colA);
            getMatAEle(matA, colA, i, arrA);

            mat *arrB = allocMatMem(rowB, 1);
            getMatBEle(matB, rowB, colB, j, arrB);

            p[(i * colB) + j] = multEle(arrA, arrB, colA);

            free(arrA);
            free(arrB);
        }
    }
}