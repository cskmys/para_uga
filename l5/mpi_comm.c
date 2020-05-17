#include<mpi.h>

#include "mpi_comm.h"


void mpi_sync_tx(double *buff, int cnt, int rxerRank){
    MPI_Ssend(buff, cnt, MPI_DOUBLE, rxerRank, 0, MPI_COMM_WORLD); // could behave as blocking call if message is too big
}

void mpi_async_tx(double *buff, int cnt, int rxerRank){
    MPI_Send(buff, cnt, MPI_DOUBLE, rxerRank, 0, MPI_COMM_WORLD); // could behave as blocking call if message is too big
}

void mpi_sync_rx(double *buff, int cnt, int txerRank){
    MPI_Recv(buff, cnt, MPI_DOUBLE, txerRank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // blocking call   
}

void mpi_chk_for_msg(int *cnt, int txerRank){
    MPI_Status status;
    MPI_Probe(txerRank, 0, MPI_COMM_WORLD, &status);
    MPI_Get_count(&status, MPI_DOUBLE, cnt);
}

