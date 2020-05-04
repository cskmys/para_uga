#include<stdio.h>
#include<assert.h>
#include "matrix_op.h"
#include "utils.h"

#define ELE_SIZ(x) sizeof(x[0])
#define NB_ELE(x) sizeof(x)/ELE_SIZ(x)

int getIdx(int nbEle, int idx){
    idx = idx % nbEle;
    if(idx < 0){
        idx = nbEle + idx;
    }
    return idx;
}

int main(int argc, char **argv){
    int n = 8;
    int p = n / 4;

    assert(n % p == 0); // important simplification we are making

    mat *A = allocMatMem(n, n);
    setIdentityMatrix(A, n);
    mat *B = allocMatMem(n, n);
    setIdentityMatrix(B, n);
    mat *C = allocMatMem(n, n);
    memsetMatrix(C, 0, n, n);
    
    for(int rank = 0; rank < p; ++rank){
        printf("For Rank: %d\n", rank);
        mat *sub_A = allocMatMem(n/p, n);
        int row_idxA = rank * n/p;
        printf("sub_A:\n");
        cpyMatrix(sub_A, A, n, n, row_idxA, n/p);
        printMatrix(sub_A, n/p, n);

        mat *sub_C = allocMatMem(n/p, n);
        memsetMatrix(sub_C, 0, n/p, n);

        for(int step = 0; step < p; ++step){
            printf("For Step: %d\n", step);
            mat *sub_B = allocMatMem(n/p, n);
            int row_idxB = getIdx(n, (rank * n/p) + (step * n/p));
            cpyMatrix(sub_B, B, n, n, row_idxB, n/p);
            printf("sub_B:\n");
            printMatrix(sub_B, n/p, n);

            mat *sub_A_A = allocMatMem(n/p, n/p);
            for (int i = 0; i < n/p; ++i){
                for (int j = 0; j < n/p; ++j){
                    int idx = getIdx( n, (((rank - step) % p) * (n/p)) + j);
                    int sub_arr_idx = (i * n) + idx;
                    sub_A_A[(i * n/p) + j] = sub_A[sub_arr_idx];
                }
            }
            printf("sub_A_A:\n");
            printMatrix(sub_A_A, n/p, n/p);
            mat *sub_C_Step = allocMatMem(n/p, n);
            size_t tmp;
            matrixMul(sub_A_A, n/p, n/p, sub_B, n/p, n, &sub_C_Step, &tmp, &tmp);
            printf("sub_C_Step:\n");
            printMatrix(sub_C_Step, n/p, n);

            sumMatrix(sub_C, sub_C_Step, n/p, n, sub_C);

            free(sub_A_A);
            free(sub_B);
        }
        printf("sub_C:\n");
        printMatrix(sub_C, n/p, n);

        size_t cCpyIdx = rank * n/p * n;
        cpyMatrix(&C[cCpyIdx], sub_C, n/p, n, 0, n/p);
        
        free(sub_C);
        free(sub_A);
    }
    printf("\n");
    
    printf("C:\n");
    printMatrix(C, n, n);

    free(A);
    free(B);
    free(C);

    return 0;
}