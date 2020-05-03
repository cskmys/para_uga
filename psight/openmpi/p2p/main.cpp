#include<iostream>
#include<stdio.h>
#include<mpi.h>
#include<math.h>

#define ELE_SIZ(x) sizeof(x[0])
#define NB_ELE(x) sizeof(x)/ELE_SIZ(x)

void mpi_printf(int rank, const char *fmt, ...){
    printf("%d : ", rank);
    va_list args;
    va_start(args, fmt);
    vprintf(fmt, args);
    va_end(args);
    return;
}

void mpi_async_tx(char *buff, int cnt, int rxRank){
    MPI_Send(buff, cnt, MPI_CHAR, rxRank, 0, MPI_COMM_WORLD);
}

void mpi_sync_tx(char *buff, int cnt, int rxRank){
    MPI_Ssend(buff, cnt, MPI_CHAR, rxRank, 0, MPI_COMM_WORLD);
}

int main(int argc, char **argv){

    MPI_Init(&argc, &argv);

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(rank == 0){
        int n = 42;
        for (int i = 1; i < size; ++i) {
            mpi_printf(rank, "Ready to tx to: %d\n", i);
            mpi_sync_tx((char*)&n, sizeof(n), i);
            mpi_printf(rank, "Data tx to: %d\n", i);
        }
    } else {
        int n;
        MPI_Recv(&n, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        mpi_printf(rank, "rxed from 0: %d\n", n);
    }


    MPI_Finalize();

    return 0;
}