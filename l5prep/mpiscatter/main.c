#include<stdio.h>
#include<stdarg.h>
#include<stdlib.h>
#include<string.h>
#include<stdbool.h>
#include<math.h>

#include<unistd.h>
#include<fcntl.h>
#include<time.h>

#include<mpi.h>

#define ELE_SIZ(x) sizeof(x[0])
#define NB_ELE(x) sizeof(x)/ELE_SIZ(x)

FILE *pFd[4] = {NULL};
MPI_Op myOp;
    
void mpi_printf(const char *fmt, ...){
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    va_list args;
    va_start(args, fmt);
    vfprintf(pFd[rank], fmt, args);
    va_end(args);
    fflush(pFd[rank]);
    return;
}

void mpi_bcast(char *buff, int cnt){
    // 0 -> txer
    MPI_Bcast(buff, cnt, MPI_CHAR, 0, MPI_COMM_WORLD);
}


void broadcast(int rank, int size){

    int value = 45;
    mpi_bcast((char*)&value, sizeof(value)); // same api with same signature to both tx and rx, pretty cool ha!
    
    if(rank == 0){
        mpi_printf("bcasted %d\n", value);   
    } else {
        mpi_printf("rxed %d from 0\n", value);    
    }
    
    MPI_Barrier(MPI_COMM_WORLD); // to understand, comment this n see

    mpi_printf("buhbye!\n");

    return;
}

void reduce(int size, int rank){

    const int itemsPerProcess = 10;
    const int count = size * itemsPerProcess;

    int *matA = (int*)malloc(count * sizeof(int));
    int *matB = (int*)malloc(count * sizeof(int));
    if(rank == 0){
        mpi_printf("matA: ");
        for (size_t i = 0; i < count; i++){
            matA[i] = i;// rand() % 10;
            mpi_printf("%d ", matA[i]);
        }
        mpi_printf("\n");

        mpi_printf("matB: ");
        for (size_t i = 0; i < count; i++){
            matB[i] = i; // rand() % 10;
            mpi_printf("%d ", matB[i]);
        }
        mpi_printf("\n");
    }

    int *localMatA = (int*)malloc(itemsPerProcess * sizeof(int));
    int *localMatB = (int*)malloc(itemsPerProcess * sizeof(int));

    // array to be scattered, nb items/proc to be scattered, data type, where to scatter the data to, how many to be recieved, data type to be recieved, scatterer, proc obj
    MPI_Scatter(matA, itemsPerProcess, MPI_INT, localMatA, itemsPerProcess, MPI_INT, 0, MPI_COMM_WORLD); // 0 -> rank of whoever's 'scattering'
    free(matA);
    MPI_Scatter(matB, itemsPerProcess, MPI_INT, localMatB, itemsPerProcess, MPI_INT, 0, MPI_COMM_WORLD); // 0 -> rank of whoever's 'scattering'
    free(matB);
    
    mpi_printf("localMatA:\n");
    for (size_t i = 0; i < itemsPerProcess; i++){
        mpi_printf("%d ", localMatA[i]);
    }
    mpi_printf("\n");

    mpi_printf("localMatB:\n");
    for (size_t i = 0; i < itemsPerProcess; i++){
        mpi_printf("%d ", localMatB[i]);
    }
    mpi_printf("\n");

    int *localMatC = (int*)malloc(itemsPerProcess * sizeof(int));
    for (size_t i = 0; i < itemsPerProcess; i++){
        localMatC[i] = localMatA[i] + localMatB[i];
    }
    free(localMatA);
    free(localMatB);
    
    mpi_printf("localMatC:\n");
    for (size_t i = 0; i < itemsPerProcess; i++){
        mpi_printf("%d ", localMatC[i]);
    }
    mpi_printf("\n");
    
    int *matC = (int*)malloc(count * sizeof(int));
    memset(matC, 0, count * sizeof(int));
    matC[0] = rank;
    matC[1] = itemsPerProcess;
    // local var to be reduced, the var to which reduced output is written, , data type, operation to be performed, scatterer, com/proc obj
    MPI_Reduce(localMatC, matC, itemsPerProcess, MPI_INT, myOp, 0, MPI_COMM_WORLD);
    free(localMatC);

    if(rank == 0){
        mpi_printf("matC:\n");
        for (size_t i = 0; i < count; i++){
            mpi_printf("%d ", matC[i]);
        }
        mpi_printf("\n");
    }
    free(matC);

    return;
}

void sumMat(int *in, int *inout, int *len, MPI_Datatype *dptr){
    int *ip = in;
    int rank = inout[0];
    int itemsPerProcess = inout[1];
    int *op = &inout[itemsPerProcess * rank];

    mpi_printf("l: %d, r: %d\n", *len, rank);
    
    for (size_t i = 0; i < *len; i++){
        mpi_printf("%d ", ip[i]);
        op[i] = ip[i];
    }
    mpi_printf("\n");
}

int main(int argc, char **argv){
    MPI_Init(&argc, &argv);

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size); // nb processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Op_create((MPI_User_function *) sumMat, false, &myOp);

    char *pipNam = strdup("./p0");
    pipNam[strlen(pipNam) - 1] = rank + '0';
    if((pFd[rank] = fopen(pipNam, "w")) < 0){
        printf("%d: fuck me!\n", rank);
        while(1);
    }

    // broadcast(size, rank);
    
    srand(time(NULL));
    reduce(size, rank);

    fclose(pFd[rank]);

    MPI_Op_free(&myOp);

    MPI_Finalize();
    return 0;
}