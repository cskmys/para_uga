#include<stdio.h>
#include<stdarg.h>
#include<stdlib.h>
#include<string.h>
#include<stdbool.h>
#include<limits.h>
#include<assert.h>

#include<unistd.h>
#include<fcntl.h>
#include<time.h>

#include<mpi.h>

#include "ring.h"
#include "matrix_op.h"
#include "utils.h"


#define ELE_SIZ(x) sizeof(x[0])
#define NB_ELE(x) sizeof(x)/ELE_SIZ(x)

#define DIM_A 16
#define COL_B 8

int getIdx(int nbEle, int idx){
    idx = idx % nbEle;
    if(idx < 0){
        idx = nbEle + idx;
    }
    return idx;
}

void reduce(int size, int rank){
    int  rA = DIM_A;
    int p = size;

    assert(rA % p == 0); // important simplification we are making
    int cA = rA;
    mat *A;

    int rB = cA;
    int cB = COL_B;
    mat *B;

    int rC = rA;
    int cC = cB;
    mat *C;
    if(rank == 0){
        A = allocMatMem(rA, cA);
        // memsetMatrix(A, 0, rA, cA);
        randsetMatrix(A, DIM_A, rA, cA);
        mpi_printf("A:\n");
        printMatrix(A, rA, cA);
        
        B = allocMatMem(rB, cB);
        randsetMatrix(B, DIM_A, rB, cB);
        mpi_printf("B:\n");
        printMatrix(B, rB, cB);

        C = allocMatMem(rC, cC);
        matrixMul(A, rA, cA, B, rB, cB, C);
        mpi_printf("C:\n");
        printMatrix(C, rC, cB);
        memsetMatrix(C, 0, rC, cC);
    }

    int rSA = rA / p;
    int cSA = cA;
    mat *sub_A = allocMatMem(rSA, cSA);
    MPI_Scatter(A, (rA / p) * cA, MPI_DOUBLE, sub_A, rSA * cSA, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if(rank == 0){
        free(A);
    }

    int rSB = rB / p;
    int cSB = cB;
    mat *sub_B = allocMatMem(rSB, cSB);
    MPI_Scatter(B, (rB / p) * cB, MPI_DOUBLE, sub_B, rSB * cSB, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if(rank == 0){
        free(B);
    }
    
    mpi_printf("sub_A:\n");
    printMatrix(sub_A, rSA, cSA);

    mpi_printf("sub_B:\n");
    printMatrix(sub_B, rSB, cSB);

    int rSC = rC / p;
    int cSC = cC;
    mat *sub_C = allocMatMem(rSC, cSC);
    memsetMatrix(sub_C, 0, rSC, cSC);

    for(int step = 0; step < p; ++step){
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
        mpi_printf("sub_A_A:\n");
        printMatrix(sub_A_A, rSAA, cSAA);

        int rSCS = rSAA;
        int cSCS = cSB;
        mat *sub_C_Step = allocMatMem(rSCS, cSCS);

        matrixMul(sub_A_A, rSAA, cSAA, sub_B, rSB, cSB, sub_C_Step);
        mpi_printf("sub_C_Step:\n");
        printMatrix(sub_C_Step, rSCS, cSCS);

        mpi_printf("cur sub_C:\n");
        printMatrix(sub_C, rSC, cSC);
        sumMatrix(sub_C, sub_C_Step, rSC, cSC, sub_C);
        mpi_printf("after sum sub_C:\n");
        printMatrix(sub_C, rSC, cSC);
        free(sub_C_Step);
        free(sub_A_A);
        
        if(rank == 0){

            ring_tx_chunk(size, rank, sub_B, rSB * cSB);
            
            int rxCnt;
            ring_rx_chunk(size, rank, sub_B, &rxCnt);
            assert(rxCnt == rSB * cSB);

            mpi_printf("Rxed:\n");
            printMatrix(sub_B, rSB, cSB);

        } else {
            mat *rxBuff = allocMatMem(rSB, cSB);
            int rxCnt;
            ring_rx_chunk(size, rank, rxBuff, &rxCnt);
            assert(rxCnt == rSB * cSB);
            
            mpi_printf("Rxed:\n");
            printMatrix(rxBuff, rSB, cSB);

            ring_tx_chunk(size, rank, sub_B, rSB * cSB);

            free(sub_B);
            sub_B = rxBuff;
        }
    }


    mpi_printf("sub_C:\n");
    printMatrix(sub_C, rSC, cSC);
    free(sub_A);
    free(sub_B);
    
    MPI_Gather(sub_C, rSC * cSC, MPI_DOUBLE, C, rSC * cSC, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    free(sub_C);

    if(rank == 0){        
        mpi_printf("matC:\n");
        printMatrix(C, rC, cC);

        free(C);
    }
    return;
}

int main(int argc, char **argv){
    MPI_Init(&argc, &argv);

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size); // nb processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    assert(!(size > MAX_NB_PROC));

    init_comm_pipe(rank);
    // broadcast(size, rank);
        
    srand(time(NULL));
    reduce(size, rank);

    deinit_comm_pipe(rank);

    MPI_Finalize();
    return 0;
}