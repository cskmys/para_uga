#include<stdio.h>
#include<assert.h>
#include<limits.h>
#include<time.h>
#include<mpi.h>

#include "utils.h"
#include "ring.h"
#include "fox.h"
#include "matrix_op.h"

typedef void(*matMulType)(int size, int rank, mat *A, int rA, int cA, mat *B, int rB, int cB, mat *C);

void evaluatePerf(int size, int rank, mat *A, int rA, int cA, mat *B, int rB, int cB, matMulType matMul, mat *CRing, mat *CSeq){
#ifdef PERF_EVAL
    for (int exp = 0; exp < NBEXPERIMENTS; exp++){
        double start;
        if(rank == 0){
            start = MPI_Wtime();
        }
        
        matMul(size, rank, A, rA, cA, B, rB, cB, CRing);
        
        if(rank == 0){
            set_exp_time(exp, MPI_Wtime() - start);
        }
    }
    if(rank == 0){
        
        double av = average_time();
        printf ("\n matrix_mul \t\t\t %.3lf seconds\n\n", av);

        for (int exp = 0 ; exp < NBEXPERIMENTS; exp++){
            double start;
                start = MPI_Wtime();
                matrixMul(A, rA, cA, B, rB, cB, CSeq);
                set_exp_time(exp, MPI_Wtime() - start);
        }
        av = average_time();
        printf ("\n seq_matrix_mul \t\t\t %.3lf seconds\n\n", av);
    }
#endif
}


void checkCorrectness(int size, int rank, mat *A, int rA, int cA, mat *B, int rB, int cB, matMulType matMul, mat *CRing, mat *CSeq){
#ifdef CHECK_CORRECTNESS

    matMul(size, rank, A, rA, cA, B, rB, cB, CRing);

    if(rank == 0){
        matrixMul(A, rA, cA, B, rB, cB, CSeq);
        int rC = rA;
        int cC = cB;
        if( cmpMatrix(CSeq, rC, cC, CRing, rC, cC) == true){
            printf("\t CORRECT matrix multiplication result \n");
        } else {
            printf("\t FAILED matrix multiplication !!! \n");
        }
    }
#endif
}

void evaluateRingPerf(int size, int rank, mat *A, int rA, int cA, mat *B, int rB, int cB, mat *CRing, mat *CSeq){
    evaluatePerf(size, rank, A, rA, cA, B, rB, cB, ring_matrix_mul, CRing, CSeq);
}

void checkRingCorrectness(int size, int rank, mat *A, int rA, int cA, mat *B, int rB, int cB, mat *CRing, mat *CSeq){
    checkCorrectness(size, rank, A, rA, cA, B, rB, cB, ring_matrix_mul, CRing, CSeq);
}

void evaluateFoxPerf(int size, int rank, mat *A, int rA, int cA, mat *B, int rB, int cB, mat *CFox, mat *CSeq){
    evaluatePerf(size, rank, A, rA, cA, B, rB, cB, fox_matrix_mul, CFox, CSeq);
}

void checkFoxCorrectness(int size, int rank, mat *A, int rA, int cA, mat *B, int rB, int cB, mat *CFox, mat *CSeq){
    checkCorrectness(size, rank, A, rA, cA, B, rB, cB, fox_matrix_mul, CFox, CSeq);
}

int main(int argc, char **argv){
    MPI_Init(&argc, &argv);

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size); // nb processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int DIM_A;
    int COL_B;
    int TYPE;
    if(argc != 4){
        if(rank == 0){
            printf("usage: %s type(0: ring, 1:fox, 2:both) matrix_A_dimension #_col_in_matrix_B\n", argv[0]);
        }
        MPI_Finalize();
        return 0;
    } else {
        TYPE = atoi(argv[1]);
        assert(TYPE >= 0);
        assert(TYPE <= 2);
        DIM_A = atoi(argv[2]);
        COL_B = atoi(argv[3]);
    }
    
#ifdef DBG_ON_FD    
    assert(!(size > MAX_NB_PROC));

    init_comm_pipe(rank);
#endif

    int  rA = DIM_A;
    int cA = rA;
    mat *A;

    int rB = cA;
    int cB = COL_B;
    mat *B;

    int rC = rA;
    int cC = cB;
    mat *CRing;
    mat *CFox;
    mat *CSeq;
    if(rank == 0){
        A = allocMatMem(rA, cA);
        randsetMatrix(A, DIM_A, rA, cA); 
        mpi_printf("A:\n");
        printMatrix(A, rA, cA);

        B = allocMatMem(rB, cB);
        randsetMatrix(B, DIM_A, rB, cB);
        mpi_printf("B:\n");
        printMatrix(B, rB, cB);

        if((TYPE == 0) || (TYPE == 2)){
            CRing = allocMatMem(rC, cC);
            memsetMatrix(CRing, 0, rC, cC);
        }
        if((TYPE == 1) || (TYPE == 2)){
            CFox = allocMatMem(rC, cC);
            memsetMatrix(CFox, 0, rC, cC);
        }

        CSeq = allocMatMem(rC, cC);
        memsetMatrix(CSeq, 0, rC, cC);
    }

    if((TYPE == 0) || (TYPE == 2)){
        evaluateRingPerf(size, rank, A, rA, cA, B, rB, cB, CRing, CSeq);
        checkRingCorrectness(size, rank, A, rA, cA, B, rB, cB, CRing, CSeq);
    }
    if((TYPE == 1) || (TYPE == 2)){
        evaluateFoxPerf(size, rank, A, rA, cA, B, rB, cB, CFox, CSeq);
        checkFoxCorrectness(size, rank, A, rA, cA, B, rB, cB, CFox, CSeq);
    }
    
    if(rank == 0){
        if((TYPE == 0) || (TYPE == 2)){
            mpi_printf("CRing:\n");
            printMatrix(CRing, rC, cC);
            free(CRing);
        }
        if((TYPE == 1) || (TYPE == 2)){
            mpi_printf("CFox:\n");
            printMatrix(CFox, rC, cC);
            free(CFox);
        }

        mpi_printf("CSeq:\n");    
        printMatrix(CSeq, rC, cC);
        free(CSeq);
    
        free(A);
        free(B);        
    }

#ifdef DBG_ON_FD
    deinit_comm_pipe(rank);
#endif

    MPI_Finalize();
    return 0;
}