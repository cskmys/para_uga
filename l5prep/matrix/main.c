#include<stdio.h>
#include<assert.h>
#include "matrix_op.h"
#include "utils.h"

#define ELE_SIZ(x) sizeof(x[0])
#define NB_ELE(x) sizeof(x)/ELE_SIZ(x)

#define NB_PROC 4
#define DIM_A 4
#define COL_B 4

int getIdx(int nbEle, int idx){
    idx = idx % nbEle;
    if(idx < 0){
        idx = nbEle + idx;
    }
    return idx;
}

int main(int argc, char **argv){
    int  rA = DIM_A;
    int p = NB_PROC;

    assert(rA % p == 0); // important simplification we are making
    int cA = rA;
    mat *A = allocMatMem(rA, cA);
    randsetMatrix(A, 10, rA, cA);
    printMatrix(A, rA, cA);
    
    int rB = cA;
    int cB = COL_B;
    mat *B = allocMatMem(rB, cB);
    randsetMatrix(B, 15, rB, cB);
    printMatrix(B, rB, cB);

    int rC = rA;
    int cC = cB;
    mat *C = allocMatMem(rC, cC);
    matrixMul(A, rA, cA, B, rB, cB, C);
    printMatrix(C, rC, cB);
    memsetMatrix(C, 0, rC, cC);

    
    for(int rank = 0; rank < p; ++rank){
        printf("For Rank: %d\n", rank);
        int rSA = rA / p;
        int cSA = cA;
        mat *sub_A = allocMatMem(rSA, cSA);
        int row_idxA = rank * rA / p;
        printf("sub_A:\n");
        cpyMatrix(sub_A, A, rA, cA, row_idxA, rSA);
        printMatrix(sub_A, rSA, cSA);

        int rSC = rC / p;
        int cSC = cC;
        mat *sub_C = allocMatMem(rSC, cSC);
        memsetMatrix(sub_C, 0, rSC, cSC);

        for(int step = 0; step < p; ++step){
            printf("For Step: %d\n", step);
            int rSB = rB / p;
            int cSB = cB; 
            mat *sub_B = allocMatMem(rSB, cSB);
            int row_idxB = getIdx(rB, (rank * rB/p) - (step * rB/p));
            cpyMatrix(sub_B, B, rB, cB, row_idxB, rSB);
            printf("sub_B:\n");
            printMatrix(sub_B, rSB, cSB);

            int rSAA = rSA;
            int cSAA = rSB;
            mat *sub_A_A = allocMatMem(rSAA, cSAA);
            for (int i = 0; i < rSAA; ++i){
                for (int j = 0; j < cSAA; ++j){
                    int idx = getIdx( cSA, (((rank - step) % p) * cSAA) + j);
                    int sub_arr_idx = (i * cSA) + idx;
                    sub_A_A[(i * rSAA) + j] = sub_A[sub_arr_idx];
                }
            }
            printf("sub_A_A:\n");
            printMatrix(sub_A_A, rSAA, cSAA);

            int rSCS = rSAA;
            int cSCS = cSB;
            mat *sub_C_Step = allocMatMem(rSCS, cSCS);

            matrixMul(sub_A_A, rSAA, cSAA, sub_B, rSB, cSB, sub_C_Step);
            printf("sub_C_Step:\n");
            printMatrix(sub_C_Step, rSCS, cSCS);

            printf("cur sub_C:\n");
            printMatrix(sub_C, rSC, cSC);
            sumMatrix(sub_C, sub_C_Step, rSC, cSC, sub_C);
            printf("after sum sub_C:\n");
            printMatrix(sub_C, rSC, cSC);

            free(sub_C_Step);
            free(sub_A_A);
            free(sub_B);
        }
        printf("sub_C:\n");
        printMatrix(sub_C, rSC, cSC);

        printf("cur C:\n");
        printMatrix(C, rC, cC);
        int cCpyIdx = rank * (rC / p) * cC;
        cpyMatrix(&C[cCpyIdx], sub_C, rSC, cSC, 0, rSC);
        printf("after cpy C:\n");
        printMatrix(C, rC, cC);

        free(sub_C);
        free(sub_A);
    }
    printf("\n");
    
    printf("C:\n");
    printMatrix(C, rC, cC);

    free(A);
    free(B);
    free(C);

    return 0;
}