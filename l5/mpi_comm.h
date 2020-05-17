#ifndef MPI_COMM_H
#define MPI_COMM_H

void mpi_sync_tx(double *buff, int cnt, int rxerRank);
void mpi_async_tx(double *buff, int cnt, int rxerRank);
void mpi_sync_rx(double *buff, int cnt, int txerRank);
void mpi_chk_for_msg(int *cnt, int txerRank);

#endif
