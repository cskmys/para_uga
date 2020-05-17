#include "fox.h"

#include<string.h>
#include<math.h>
#include<mpi.h>

#include "cart.h"
#include "mpi_comm.h"

#include "matrix_op.h"

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
void orgMatrix(mat *M, int rM, int cM, mat *mM, int sizSqrtI){
    int size = sizSqrtI * sizSqrtI;
    int rMM = rM;
    int cMM = cM;
    for (size_t i = 0; i < size; i++){
        int rtM = rM / sizSqrtI;
        int ctM = cM / sizSqrtI;
        mat *tM = allocMatMem(rtM, ctM);
        int idx = (int)((rMM * cMM)/size) * i;
        buildMatrix(M, rM, cM, tM, i, sizSqrtI);
        memcpy(&mM[ idx ], tM, rtM * ctM * sizeof(mat));
        free(tM);
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
void reorgMatrix(mat *mM, int rMM, int cMM, mat *M, int sizSqrtI){
    int size = sizSqrtI * sizSqrtI;
    int rM = rMM;
    int cM = cMM;
    for (size_t i = 0; i < size; i++){
        int rtM = rM / sizSqrtI;
        int ctM = cM / sizSqrtI;
        mat *tM = allocMatMem(rtM, ctM);
        int idx = (int)((rMM * cMM)/size) * i;
        memcpy(tM, &mM[ idx ], rtM * ctM * sizeof(mat));
        rebuildMatrix(tM, rtM, ctM, M, i, sizSqrtI);
        free(tM);
    }
}

void fox_matrix_mul(int size, int rank, mat *A, int rA, int cA, mat *B, int rB, int cB, mat *C){
    int sizSqrtI = (int)sqrt(size);
    double sizSqrtD = sqrt(size);
    assert((double)sizSqrtI == sizSqrtD); // we need a perfect square

    assert(rA == cA); // pre-condition
    assert(cA == rB); // validaity of operation
    assert(rA % sizSqrtI == 0); // important simplification we expect
    assert(cA % sizSqrtI == 0);
    assert(cB % sizSqrtI == 0);
    
    int rMA = rA;
    int cMA = cA;
    mat *mA = allocMatMem(rMA, cMA);
    int rMB = rB;
    int cMB = cB;
    mat *mB = allocMatMem(rMB, cMB);
    int rMC = rMA;
    int cMC = cMB;
    mat *mC = allocMatMem(rMC, cMC);
    if(rank == 0){
        orgMatrix(A, rA, cA, mA, sizSqrtI);
        orgMatrix(B, rB, cB, mB, sizSqrtI);
    }

    int rC = rA;
    int cC = cB;

    int rSA = rMA / sizSqrtI;
    int cSA = cMA / sizSqrtI;
    mat *sub_A = allocMatMem(rSA, cSA);
    MPI_Scatter(mA, (rMA / sizSqrtI) * (cMA / sizSqrtI), MPI_MAT, sub_A, rSA * cSA, MPI_MAT, 0, MPI_COMM_WORLD);
    if(rank == 0){
        free(mA);
    }

    int rSB = rB / sizSqrtI;
    int cSB = cB / sizSqrtI;
    mat *sub_B = allocMatMem(rSB, cSB);
    MPI_Scatter(mB, (rMB / sizSqrtI) * (cMB / sizSqrtI), MPI_MAT, sub_B, rSB * cSB, MPI_MAT, 0, MPI_COMM_WORLD);
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
    
    MPI_Comm comm1dRow[sizSqrtI];
    splitCartByRow(comm2d, &comm1dRow[coord2d[0]]);
    int coord1dRow[NB_DIM];
    getRowCoord(comm1dRow[coord2d[0]],  coord1dRow);

    MPI_Comm comm1dCol[sizSqrtI];
    splitCartByCol(comm2d, &comm1dCol[coord2d[1]]);
    int coord1dCol[NB_DIM];
    getColCoord(comm1dCol[coord2d[1]],  coord1dCol);

    for (size_t rStep = 0; rStep < sizSqrtI; rStep++){
        int rSAA = rSA;
        int cSAA = cSA;
        mat *sub_A_A = allocMatMem(rSAA, cSAA);
        memsetMatrix(sub_A_A, 0.0, rSAA, cSAA);
        int rIdx = (rStep + (int)(rank/sizSqrtI)) % sizSqrtI;
        if(getRowRank(comm1dRow[coord2d[0]]) == rIdx){
            memcpy(sub_A_A, sub_A, rSAA * cSAA * sizeof(mat));
        }
        MPI_Bcast(sub_A_A, rSAA * cSAA, MPI_MAT, rIdx, comm1dRow[coord2d[0]]);
        
        int rSBB = rSB;
        int cSBB = cSB;
        mat *sub_B_B = allocMatMem(rSBB, cSBB);
        memsetMatrix(sub_B_B, 0.0, rSBB, cSBB);
        int txIdx = (getColRank(comm1dCol[coord2d[1]]) + sizSqrtI - rStep) % sizSqrtI;
        int rxIdx = (getColRank(comm1dCol[coord2d[1]]) + rStep ) % sizSqrtI;
        
        MPI_Send(sub_B, rSB * cSB, MPI_MAT, txIdx, 0, comm1dCol[coord2d[1]]);
        MPI_Status status;
        MPI_Recv(sub_B_B, rSBB * cSBB, MPI_MAT, rxIdx, 0, comm1dCol[coord2d[1]], &status);
        
        int rSCS = rSA;
        int cSCS = cSB;
        mat *sub_C_Step = allocMatMem(rSCS, cSCS);

        matrixMul(sub_A_A, rSAA, cSAA, sub_B_B, rSBB, cSBB, sub_C_Step);

        free(sub_A_A);
        free(sub_B_B);

        sumMatrix(sub_C, sub_C_Step, rSC, cSC, sub_C);
        free(sub_C_Step);
    }

    delCol(&comm1dCol[coord2d[1]]);
    delRow(&comm1dRow[coord2d[0]]);    
    delCartTopology(&comm2d);

    free(sub_B);
    free(sub_A);

    MPI_Gather(sub_C, rSC * cSC, MPI_MAT, mC, rSC * cSC, MPI_MAT, 0, MPI_COMM_WORLD);
    if(rank == 0){
        reorgMatrix(mC, rMC, cMC, C, sizSqrtI);
    }
    free(sub_C);
}