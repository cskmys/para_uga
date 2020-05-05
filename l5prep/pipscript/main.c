#include<stdio.h>
#include<stdarg.h>
#include<stdlib.h>
#include<string.h>
#include<stdbool.h>
#include<math.h>
#include<assert.h>

#include<unistd.h>
#include<fcntl.h>
#include<time.h>

#include<mpi.h>

#include "ring.h"

#define ELE_SIZ(x) sizeof(x[0])
#define NB_ELE(x) sizeof(x)/ELE_SIZ(x)


void reduce(int size, int rank){

    const int itemsPerProcess = 10;
    const int count = size * itemsPerProcess;

    double *matA;
    double *matB;
    if(rank == 0){
        matA = (double*)malloc(count * sizeof(double));
        matB = (double*)malloc(count * sizeof(double));
        mpi_printf("matA: ");
        for (size_t i = 0; i < count; i++){
            matA[i] = i;// rand() % 10;
            mpi_printf("%f ", matA[i]);
        }
        mpi_printf("\n");

        mpi_printf("matB: ");
        for (size_t i = 0; i < count; i++){
            matB[i] = i; // rand() % 10;
            mpi_printf("%f ", matB[i]);
        }
        mpi_printf("\n");
    }

    double *localMatA = (double*)malloc(itemsPerProcess * sizeof(double));
    double *localMatB = (double*)malloc(itemsPerProcess * sizeof(double));

    // array to be scattered, nb items/proc to be scattered, data type, where to scatter the data to, how many to be recieved, data type to be recieved, scatterer, proc obj
    MPI_Scatter(matA, itemsPerProcess, MPI_DOUBLE, localMatA, itemsPerProcess, MPI_DOUBLE, 0, MPI_COMM_WORLD); // 0 -> rank of whoever's 'scattering'
    if(rank == 0){
        free(matA);
    }
    MPI_Scatter(matB, itemsPerProcess, MPI_DOUBLE, localMatB, itemsPerProcess, MPI_DOUBLE, 0, MPI_COMM_WORLD); // 0 -> rank of whoever's 'scattering'
    if(rank == 0){
        free(matB);
    }
    
    mpi_printf("localMatA:\n");
    for (size_t i = 0; i < itemsPerProcess; i++){
        mpi_printf("%f ", localMatA[i]);
    }
    mpi_printf("\n");

    mpi_printf("localMatB:\n");
    for (size_t i = 0; i < itemsPerProcess; i++){
        mpi_printf("%f ", localMatB[i]);
    }
    mpi_printf("\n");


    if(rank == 0){

        ring_tx_chunk(size, rank, localMatB, itemsPerProcess);
        
        int rxCnt;
        ring_rx_chunk(size, rank, localMatB, &rxCnt);
        assert(rxCnt == itemsPerProcess);

        mpi_printf("Rxed:\n");
        for (size_t i = 0; i < itemsPerProcess; i++){
            mpi_printf("%f ", localMatB[i]);
        }
        mpi_printf("\n");

    } else {
        double *rxBuff = (double*)malloc(itemsPerProcess * sizeof(double));
        int rxCnt;
        ring_rx_chunk(size, rank, rxBuff, &rxCnt);
        assert(rxCnt == itemsPerProcess);
        
        mpi_printf("Rxed:\n");
        for (size_t i = 0; i < itemsPerProcess; i++){
            mpi_printf("%f ", rxBuff[i]);
        }
        mpi_printf("\n");

        ring_tx_chunk(size, rank, localMatB, itemsPerProcess);

        memcpy(localMatB, rxBuff, itemsPerProcess * sizeof(double));
    }

    double *localMatC = (double*)malloc(itemsPerProcess * sizeof(double));
    for (size_t i = 0; i < itemsPerProcess; i++){
        localMatC[i] = localMatA[i] + localMatB[i];
    }
    free(localMatA);
    free(localMatB);
    
    mpi_printf("localMatC:\n");
    for (size_t i = 0; i < itemsPerProcess; i++){
        mpi_printf("%f ", localMatC[i]);
    }
    mpi_printf("\n");
    
    double *matC;
    if(rank == 0){
        matC = (double*)malloc(count * sizeof(double));
        memset(matC, 0, count * sizeof(double));
    }
    
    MPI_Gather(localMatC, itemsPerProcess, MPI_DOUBLE, matC, itemsPerProcess, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if(rank == 0){
        free(localMatC);

        mpi_printf("matC:\n");
        for (size_t i = 0; i < count; i++){
            mpi_printf("%f ", matC[i]);
        }
        mpi_printf("\n");

        free(matC);
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