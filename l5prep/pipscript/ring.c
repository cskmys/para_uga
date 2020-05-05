#include<stdio.h>
#include<stdarg.h>
#include<string.h>
#include<stdlib.h>

#include<mpi.h>

#include "ring.h"

FILE *pFd[MAX_NB_PROC] = {NULL};

void init_comm_pipe(int rank){
    char *filNam = strdup("f0");
    filNam[strlen(filNam) - 1] = rank + '0';

    if((pFd[rank] = fopen(filNam, "wb")) < 0){
        printf("%d: open fuck me!\n", rank);
        free(filNam);
        while(1);
    }
    free(filNam);
}

void deinit_comm_pipe(int rank){
    fclose(pFd[rank]);
}

void mpi_fprintf(FILE *f, const char *fmt, va_list args){
    
    vfprintf(f, fmt, args);
    
    fflush(f);
    return;
}

void mpi_printf(const char *fmt, ...){
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    va_list args;
    va_start(args, fmt);

    mpi_fprintf(pFd[rank], fmt, args);

    va_end(args);
    return;
}

void mpi_fprintf_rank(FILE *f, const char *fmt, ...){
    va_list args;
    va_start(args, fmt);
    
    mpi_fprintf(f, fmt, args);
    
    va_end(args);
    return;
}

void mpi_scr_printf(const char *fmt, ...){
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    mpi_fprintf_rank(stdout, "%d : ", rank);
    
    va_list args;
    va_start(args, fmt);
    
    mpi_fprintf(stdout, fmt, args);
    
    va_end(args);
    return;
}

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

void ring_tx(int myRank, int nbProc, double *buff, int cnt){
    mpi_async_tx(buff, cnt, (myRank + 1) % nbProc);
}

void ring_stx(int myRank, int nbProc, double *buff, int cnt){
    mpi_sync_tx(buff, cnt, (myRank + 1) % nbProc);
}

void ring_rx(int myRank, int nbProc, double *buff, int cnt){
    mpi_sync_rx(buff, cnt, (myRank - 1) < 0 ? (nbProc - 1) : (myRank - 1));
}

void ring_chk_for_msg(int myRank, int nbProc, int *cnt){
    mpi_chk_for_msg(cnt, (myRank - 1) < 0 ? (nbProc - 1) : (myRank - 1));
}

void ring_rx_chunk(int nbProc, int myRank, double *buff, int *nbArrEle){
    int nbEle;
    ring_chk_for_msg(myRank, nbProc, &nbEle);
    ring_rx(myRank, nbProc, buff, nbEle);
    *nbArrEle = nbEle;
}

void ring_tx_chunk(int nbProc, int myRank, double *buff, int nbArrEle){
    ring_tx(myRank, nbProc, buff, nbArrEle);
}

// void broadcast_ring(int myRank, int nbArrEle, int nbProc){
//     if(myRank == 0){
//         // send 1st to everyone
//         // then do my computation
//         // wait for result from everyone
//         int arr[nbArrEle]; // simple case of nbProc = array size
//         for (size_t i = 0; i < NB_ELE(arr); i++){
//             arr[i] = i;
//         }
//         mpi_printf(myRank, "inited\n");
//         for(int gi = 1; gi < nbProc; ++gi){ // cause 0 is txer
//             ring_stx(myRank, nbProc, (char*)&arr[gi * NB_ELE(arr)/nbProc], NB_ELE(arr)/nbProc * sizeof(int), gi);
//         }
//         mpi_printf(myRank, "txed\n");
//         proc_my_chunk(nbArrEle, nbProc, arr);
//         mpi_printf(myRank, "proced\n");

//         for(int gi = 1; gi < nbProc; ++gi){ // cause 0 is rxer
//             int buff[NB_ELE(arr)/nbProc];
//             int tag;
//             mpi_printf(myRank, "wait\n");
//             ring_rx_chunk(nbArrEle, nbProc, myRank, buff, &tag);
//             assert(tag == 0);
//             mpi_printf(myRank, "rxed1\n");
//             int arrIdx = tag * NB_ELE(arr)/nbProc;
//             memcpy((char*)&arr[arrIdx], buff, sizeof(buff));
//         }
//         mpi_printf(myRank, "Array output:");
//         for (size_t i = 0; i < NB_ELE(arr); i++){
//             printf("%d ", arr[i]);
//         }
//         printf("\n");
//     } else {
//         // recieve 1st, 
//         // if your's do computation and then send
//         // else just send
//         while(1){
//             int buff[nbArrEle/nbProc];
//             int tag;
//             mpi_printf(myRank, "wait\n");
//             ring_rx_chunk(nbArrEle, nbProc, myRank, buff, &tag);
//             mpi_printf(myRank, "rxed1\n");

//             if(tag == myRank){
//                 proc_my_chunk(nbArrEle, nbProc, buff);
                
//                 ring_tx(myRank, nbProc, (char*)buff, sizeof(buff), 0);
//                 mpi_printf(myRank, "protxed\n");
//                 continue;
//             }

//             ring_tx(myRank, nbProc, (char*)buff, sizeof(buff), tag);
//             mpi_printf(myRank, "txed\n");
//         }
//     }
// }