#include "fox.h"

#include<string.h>
#include<math.h>
#include<mpi.h>

#include "cart.h"
#include "mpi_comm.h"

#include "matrix_op.h"

// void fox_matrix_mul(int size, int rank, mat *A, int rA, int cA, mat *B, int rB, int cB, mat *C){
//     assert(rA == cA); // pre-condition
//     assert(rA % size == 0); // important simplification we expect
//     assert(rB == cA); // validity of operation

//     int p = size;

//     int rC = rA;
//     int cC = cB;
    
//     int rSA = rA / p;
//     int cSA = cA;
//     mat *sub_A = allocMatMem(rSA, cSA);
//     MPI_Scatter(A, (rA / p) * cA, MPI_DOUBLE, sub_A, rSA * cSA, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
//     int rSB = rB / p;
//     int cSB = cB;
//     mat *sub_B = allocMatMem(rSB, cSB);
//     MPI_Scatter(B, (rB / p) * cB, MPI_DOUBLE, sub_B, rSB * cSB, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
//     int rSC = rC / p;
//     int cSC = cC;
//     mat *sub_C = allocMatMem(rSC, cSC);
//     memsetMatrix(sub_C, 0, rSC, cSC);

//     for(int step = 0; step < p; ++step){
//         int rSAA = rSA;
//         int cSAA = rSB;
//         mat *sub_A_A = allocMatMem(rSAA, cSAA);
//         for (int i = 0; i < rSAA; ++i){
//             for (int j = 0; j < cSAA; ++j){
//                 int idx = getIdx( cSA, (((rank - step) % p) * cSAA) + j);
//                 int sub_arr_idx = (i * cSA) + idx;
//                 sub_A_A[(i * rSAA) + j] = sub_A[sub_arr_idx];
//             }
//         }

//         int rSCS = rSAA;
//         int cSCS = cSB;
//         mat *sub_C_Step = allocMatMem(rSCS, cSCS);

//         matrixMul(sub_A_A, rSAA, cSAA, sub_B, rSB, cSB, sub_C_Step);

//         sumMatrix(sub_C, sub_C_Step, rSC, cSC, sub_C);
//         free(sub_C_Step);
//         free(sub_A_A);
        
//         if(rank == 0){
//             bcast_tx_chunk(size, rank, sub_B, rSB * cSB);
            
//             int rxCnt;
//             bcast_rx_chunk(size, rank, sub_B, &rxCnt);
//             assert(rxCnt == rSB * cSB);
//         } else {
//             mat *rxBuff = allocMatMem(rSB, cSB);
//             int rxCnt;
//             bcast_rx_chunk(size, rank, rxBuff, &rxCnt);
//             assert(rxCnt == rSB * cSB);
            
//             bcast_tx_chunk(size, rank, sub_B, rSB * cSB);

//             free(sub_B);
//             sub_B = rxBuff;
//         }
//     }


//     free(sub_A);
//     free(sub_B);
    
//     MPI_Gather(sub_C, rSC * cSC, MPI_DOUBLE, C, rSC * cSC, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//     free(sub_C);

//     return;
// }

// void fox_mat_mul_prep(int size, int rank, mat *A, int rA, int cA, mat *B, int rB, int cB, mat *C){
//     int sizSqrtI = (int)sqrt(size);
//     double sizSqrtD = sqrt(size);
//     assert((double)sizSqrtI == sizSqrtD); // we need a perfect square

//     assert(rA == cA); // pre-condition
//     assert(rA % sizSqrtI == 0); // important simplification we expect
//     assert(rB == cA); // validity of operation
    
//     UNUSED(A);
//     UNUSED(B);
//     UNUSED(C);
//     UNUSED(cB);

//     MPI_Comm comm2d;
//     createCartTopology(sizSqrtI, sizSqrtI, &comm2d);
//     int coord2d[NB_DIM];
//     getCartCoord(comm2d,  coord2d);
    
//     mpi_printf("coord within cart: x:%d, y:%d\n", coord2d[0], coord2d[1]);
//     mpi_printf("rank within cart: %d\n", getCartRank(comm2d));

//     MPI_Comm comm1dRow;
//     splitCartByRow(comm2d, &comm1dRow);
//     int coord1d[NB_DIM];
//     getRowCoord(comm1dRow,  coord1d);

//     mpi_printf("coord within row: x:%d, y:%d\n", coord1d[0], coord1d[1]);
//     mpi_printf("rank within row: %d\n", getRowRank(comm1dRow));

//     MPI_Comm comm1dCol;
//     splitCartByCol(comm2d, &comm1dCol);
//     getColCoord(comm1dCol,  coord1d);

//     mpi_printf("coord within col: x:%d, y:%d\n", coord1d[0], coord1d[1]);
//     mpi_printf("rank within col: %d\n", getColRank(comm1dCol));

//     delCol(&comm1dCol);
//     delRow(&comm1dRow);
//     delCartTopology(&comm2d);
// }
/*
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
*/
void buildMatrix(mat *m, int rM, int cM, mat *sub_M, int rank, int sizSqrtI){
    int rSM = rM/sizSqrtI;
    int cSM = cM/sizSqrtI;
    for (int i = 0; i < rSM; ++i){
        for (int j = 0; j < cSM; ++j){
            int idx = getIdx( cM, ((rank % (sizSqrtI * sizSqrtI)) * cSM) + j);
            int arr_idx =  ( (i + ( ((int)(rank/sizSqrtI)) * ((int)(rM/sizSqrtI))) ) * cM) + idx;
            sub_M[(i * cSM) + j] = m[arr_idx];
        }
    }
}

void rebuildMatrix(mat *sub_M, int rSM, int cSM, mat *m, int rank, int sizSqrtI){
    int rM = rSM * sizSqrtI;
    int cM = cSM * sizSqrtI;
    for (int i = 0; i < rSM; ++i){
        for (int j = 0; j < cSM; ++j){
            int idx = getIdx( cM, ((rank % (sizSqrtI * sizSqrtI)) * cSM) + j);
            int arr_idx =  ( (i + ( ((int)(rank/sizSqrtI)) * ((int)(rM/sizSqrtI))) ) * cM) + idx;
            m[arr_idx] = sub_M[(i * cSM) + j];
        }
    }
}

void fox_mat_mul_prep(int size, int rank, mat *A, int rA, int cA, mat *B, int rB, int cB, mat *C){
    int sizSqrtI = (int)sqrt(size);
    double sizSqrtD = sqrt(size);
    assert((double)sizSqrtI == sizSqrtD); // we need a perfect square

    assert(rA == cA); // pre-condition
    assert(rA % sizSqrtI == 0); // important simplification we expect
    assert(rB == cA); // validity of operation
    
    int rMA = rA;
    int cMA = cA;
    mat *mA = allocMatMem(rMA, cMA);
    int rMB = rB;
    int cMB = cB;
    mat *mB = allocMatMem(rMB, cMB);
    int rMC = rA;
    int cMC = cB;
    mat *mC = allocMatMem(rMC, cMC);
    if(rank == 0){
        for (size_t i = 0; i < size; i++){
            int rtA = rA / sizSqrtI;
            int ctA = cA / sizSqrtI;
            mat *tA = allocMatMem(rtA, ctA);
            buildMatrix(A, rA, cA, tA, i, sizSqrtI);
            memcpy(&mA[ ((int)(rMA/size)) * i * cMA ], tA, rtA * ctA * sizeof(mat));
            free(tA);
        }
        for (size_t i = 0; i < size; i++){
            int rtB = rB / sizSqrtI;
            int ctB = cB / sizSqrtI;
            mat *tB = allocMatMem(rtB, ctB);
            buildMatrix(B, rB, cB, tB, i, sizSqrtI);
            memcpy(&mB[ ((int)(rMB/size)) * i * cMB ], tB, rtB * ctB * sizeof(mat));
            free(tB);
        }
    }

    int rC = rA;
    int cC = cB;

    int rSA = rA / sizSqrtI;
    int cSA = cA / sizSqrtI;
    mat *sub_A = allocMatMem(rSA, cSA);
    MPI_Scatter(mA, rMA * cMA / size, MPI_DOUBLE, sub_A, rSA * cSA, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    printMatrix(sub_A, rSA, cSA);
    if(rank == 0){
        free(mA);
    }

    int rSB = rB / sizSqrtI;
    int cSB = cB / sizSqrtI;
    mat *sub_B = allocMatMem(rSB, cSB);
    MPI_Scatter(mB, rMB * cMB / size, MPI_DOUBLE, sub_B, rSB * cSB, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    printMatrix(sub_B, rSB, cSB);
    if(rank == 0){
        free(mB);
    }

    int rSC = rC / sizSqrtI;
    int cSC = cC / sizSqrtI;
    mat *sub_C = allocMatMem(rSC, cSC);
    memsetMatrix(sub_C, 0, rSC, cSC);


    MPI_Comm comm2d;
    createCartTopology(sizSqrtI, sizSqrtI, &comm2d);
    int coord2d[NB_DIM];
    getCartCoord(comm2d,  coord2d);
    
    mpi_printf("coord within cart: x:%d, y:%d\n", coord2d[0], coord2d[1]);
    mpi_printf("rank within cart: %d\n", getCartRank(comm2d));

    MPI_Comm comm1dRow[sizSqrtI];
    splitCartByRow(comm2d, &comm1dRow[coord2d[0]]);
    int coord1dRow[NB_DIM];
    getRowCoord(comm1dRow[coord2d[0]],  coord1dRow);

    mpi_printf("coord within row: %d\n", coord1dRow[0]);
    mpi_printf("rank within row: %d\n", getRowRank(comm1dRow[coord2d[0]]));    
    
    MPI_Comm comm1dCol[sizSqrtI];
    splitCartByCol(comm2d, &comm1dCol[coord2d[1]]);
    int coord1dCol[NB_DIM];
    getColCoord(comm1dCol[coord2d[1]],  coord1dCol);

    mpi_printf("coord within col: %d\n", coord1dCol[0]);
    mpi_printf("rank within col: %d\n", getColRank(comm1dCol[coord2d[1]]));

    int rowBcastVal[NB_ELE(coord2d)];
    memcpy(rowBcastVal, coord2d, sizeof(rowBcastVal));
    MPI_Bcast(rowBcastVal, 2, MPI_INT, 0, comm1dRow[coord2d[0]]);
    mpi_printf("row bcast: %d %d\n", rowBcastVal[0], rowBcastVal[1]);

    int colBcastVal[NB_ELE(coord2d)];
    memcpy(colBcastVal, coord2d, sizeof(colBcastVal));
    MPI_Bcast(colBcastVal, 2, MPI_INT, 0, comm1dCol[coord2d[1]]);
    mpi_printf("col bcast: %d %d\n", colBcastVal[0], colBcastVal[1]);

    delCol(&comm1dCol[coord2d[1]]);
    delRow(&comm1dRow[coord2d[0]]);    
    delCartTopology(&comm2d);

    free(sub_B);
    free(sub_A);

    MPI_Gather(sub_C, rSC * cSC, MPI_DOUBLE, mC, rSC * cSC, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if(rank == 0){
        for (size_t i = 0; i < size; i++){
            int rtC = rC / sizSqrtI;
            int ctC = cC / sizSqrtI;
            mat *tC = allocMatMem(rtC, ctC);
            memcpy(tC, &mC[ ((int)(rMC/size)) * i * cMC ], rtC * ctC * sizeof(mat));
            rebuildMatrix(tC, rtC, ctC, C, i, sizSqrtI);
            free(tC);
        }
    }
    free(sub_C);
}