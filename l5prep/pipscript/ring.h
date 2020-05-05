#ifndef RING_H
#define RING_H

#define MAX_NB_PROC 64

#define ELE_SIZ(x) sizeof(x[0])
#define NB_ELE(x) sizeof(x)/ELE_SIZ(x)

void init_comm_pipe(int rank);
void deinit_comm_pipe(int rank);

void mpi_printf(const char *fmt, ...);
void mpi_scr_printf(const char *fmt, ...);

void ring_rx_chunk(int nbProc, int myRank, double *buff, int *nbArrEle);
void ring_tx_chunk(int nbProc, int myRank, double *buff, int nbArrEle);

#endif