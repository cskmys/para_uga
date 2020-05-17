#include<stdio.h>
#include<assert.h>
#include<limits.h>
#include<time.h>
#include<mpi.h>

#include "utils.h"
#include "ring.h"
#include "fox.h"
#include "matrix_op.h"


void evaluatePerf(int size, int rank, mat *A, int rA, int cA, mat *B, int rB, int cB, mat *CRing, mat *CSeq){
#ifdef PERF_EVAL
    for (int exp = 0; exp < NBEXPERIMENTS; exp++){
        double start;
        if(rank == 0){
            start = MPI_Wtime();
        }
        
        ring_matrix_mul(size, rank, A, rA, cA, B, rB, cB, CRing);
        
        if(rank == 0){
            set_exp_time(exp, MPI_Wtime() - start);
        }
    }
    if(rank == 0){
        
        double av = average_time();
        printf ("\n ring_matrix_mul \t\t\t %.3lf seconds\n\n", av);

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


void checkCorrectness(int size, int rank, mat *A, int rA, int cA, mat *B, int rB, int cB, mat *CRing, mat *CSeq){
#ifdef CHECK_CORRECTNESS
        
    ring_matrix_mul(size, rank, A, rA, cA, B, rB, cB, CRing);

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

int main(int argc, char **argv){
    MPI_Init(&argc, &argv);

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size); // nb processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int DIM_A;
    int COL_B;
    if(argc != 3){
        if(rank == 0){
            printf("usage: %s matrix_A_dimension #_col_in_matrix_B\n", argv[0]);
        }
        MPI_Finalize();
        return 0;
    } else {
        DIM_A = atoi(argv[1]);
        COL_B = atoi(argv[2]);
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

        CRing = allocMatMem(rC, cC);
        memsetMatrix(CRing, 0, rC, cC);

        CFox = allocMatMem(rC, cC);
        memsetMatrix(CFox, 0, rC, cC);

        CSeq = allocMatMem(rC, cC);
        memsetMatrix(CSeq, 0, rC, cC);
    }

    fox_mat_mul_prep(size, rank, A, rA, cA, B, rB, cB, CFox);
    // evaluatePerf(size, rank, A, rA, cA, B, rB, cB, CRing, CSeq);
    // checkCorrectness(size, rank, A, rA, cA, B, rB, cB, CRing, CSeq);
    
    if(rank == 0){
        mpi_printf("CRing:\n");
        printMatrix(CRing, rC, cC);

        mpi_printf("CFox:\n");
        printMatrix(CFox, rC, cC);

        mpi_printf("CSeq:\n");    
        printMatrix(CSeq, rC, cC);
    
        free(A);
        free(B);
        free(CSeq);
        free(CRing);
    }

#ifdef DBG_ON_FD
    deinit_comm_pipe(rank);
#endif

    MPI_Finalize();
    return 0;
}