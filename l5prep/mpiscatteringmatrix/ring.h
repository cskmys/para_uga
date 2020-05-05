#ifndef RING_H
#define RING_H

void ring_rx_chunk(int nbProc, int myRank, double *buff, int *nbArrEle);
void ring_tx_chunk(int nbProc, int myRank, double *buff, int nbArrEle);

#endif