#include<mpi.h>

#include "fox.h"
#include "mpi_comm.h"

#include "matrix_op.h"

void fox_matrix_mul(int size, int rank, mat *A, int rA, int cA, mat *B, int rB, int cB, mat *C){
    assert(rA == cA); // pre-condition
    assert(rA % size == 0); // important simplification we expect
    assert(rB == cA); // validity of operation

    int p = size;

    int rC = rA;
    int cC = cB;
    
    int rSA = rA / p;
    int cSA = cA;
    mat *sub_A = allocMatMem(rSA, cSA);
    MPI_Scatter(A, (rA / p) * cA, MPI_DOUBLE, sub_A, rSA * cSA, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    int rSB = rB / p;
    int cSB = cB;
    mat *sub_B = allocMatMem(rSB, cSB);
    MPI_Scatter(B, (rB / p) * cB, MPI_DOUBLE, sub_B, rSB * cSB, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
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

        int rSCS = rSAA;
        int cSCS = cSB;
        mat *sub_C_Step = allocMatMem(rSCS, cSCS);

        matrixMul(sub_A_A, rSAA, cSAA, sub_B, rSB, cSB, sub_C_Step);

        sumMatrix(sub_C, sub_C_Step, rSC, cSC, sub_C);
        free(sub_C_Step);
        free(sub_A_A);
        
        if(rank == 0){
            bcast_tx_chunk(size, rank, sub_B, rSB * cSB);
            
            int rxCnt;
            bcast_rx_chunk(size, rank, sub_B, &rxCnt);
            assert(rxCnt == rSB * cSB);
        } else {
            mat *rxBuff = allocMatMem(rSB, cSB);
            int rxCnt;
            bcast_rx_chunk(size, rank, rxBuff, &rxCnt);
            assert(rxCnt == rSB * cSB);
            
            bcast_tx_chunk(size, rank, sub_B, rSB * cSB);

            free(sub_B);
            sub_B = rxBuff;
        }
    }


    free(sub_A);
    free(sub_B);
    
    MPI_Gather(sub_C, rSC * cSC, MPI_DOUBLE, C, rSC * cSC, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    free(sub_C);

    return;
}