#include<iostream>
#include<stdio.h>
#include<mpi.h>
#include<math.h>
#include<assert.h>

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

void mpi_async_tx(char *buff, int cnt, int tag, int rxerRank){
    MPI_Send(buff, cnt, MPI_CHAR, rxerRank, tag, MPI_COMM_WORLD); // could behave as blocking call if message is too big
}

void mpi_sync_rx(char *buff, int cnt, int tag, int txerRank){
    MPI_Recv(buff, cnt, MPI_CHAR, txerRank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // blocking call   
}

void mpi_chk_for_msg(int *cnt, int *tag, int txerRank){
    mpi_printf(txerRank + 1, "prob\n");
    MPI_Status status;
    MPI_Probe(txerRank, 0, MPI_COMM_WORLD, &status);
    MPI_Get_count(&status, MPI_CHAR, cnt);
    *tag = status.MPI_TAG;
}

void ring_tx(int myRank, int nbProc, char *buff, int cnt, int tag){
    mpi_async_tx(buff, cnt, tag, (myRank + 1) % nbProc);
}

void ring_rx(int myRank, int nbProc, char *buff, int cnt, int tag){
    mpi_sync_rx(buff, cnt, tag, (myRank - 1) < 0 ? 0 : (myRank - 1));
}

void ring_chk_for_msg(int myRank, int nbProc, int *cnt, int *tag){
    mpi_chk_for_msg(cnt, tag, (myRank - 1) < 0 ? 0 : (myRank - 1));
}

void ring_rx_chunk(int nbArrEle, int nbProc, int myRank, int *buff, int *tag){
    int cnt;
    ring_chk_for_msg(myRank, nbProc, &cnt, tag);
    assert(cnt == nbArrEle/nbProc * sizeof(int));
    ring_rx(myRank, nbProc, (char*)buff, cnt, *tag);
}

int unit_oper(int val){
    return (val * 2);
}

void proc_my_chunk(int nbArrEle, int nbProc, int *buff){
    for(int li = 0; li < nbArrEle/nbProc; ++li){
        buff[li] = unit_oper(buff[li]);
    }
}

void broadcast_ring(int myRank, int nbArrEle, int nbProc){
    if(myRank == 0){
        // send 1st to everyone
        // then do my computation
        // wait for result from everyone
        int arr[nbArrEle]; // simple case of nbProc = array size
        for (size_t i = 0; i < NB_ELE(arr); i++){
            arr[i] = i;
        }
        mpi_printf(myRank, "inited\n");
        for(int gi = 1; gi < nbProc; ++gi){ // cause 0 is txer
            ring_tx(myRank, nbProc, (char*)&arr[gi * NB_ELE(arr)/nbProc], NB_ELE(arr)/nbProc * sizeof(int), gi);
        }
        mpi_printf(myRank, "txed\n");
        proc_my_chunk(nbArrEle, nbProc, arr);
        mpi_printf(myRank, "proced\n");

        for(int gi = 1; gi < nbProc; ++gi){ // cause 0 is rxer
            int buff[NB_ELE(arr)/nbProc];
            int tag;
            mpi_printf(myRank, "wait\n");
            ring_rx_chunk(nbArrEle, nbProc, myRank, buff, &tag);
            assert(tag == 0);
            mpi_printf(myRank, "rxed1\n");
            int arrIdx = tag * NB_ELE(arr)/nbProc;
            memcpy((char*)&arr[arrIdx], buff, sizeof(buff));
        }
        mpi_printf(myRank, "Array output:");
        for (size_t i = 0; i < NB_ELE(arr); i++){
            printf("%d ", arr[i]);
        }
        printf("\n");
    } else {
        // recieve 1st, 
        // if your's do computation and then send
        // else just send
        while(1){
            int buff[nbArrEle/nbProc];
            int tag;
            mpi_printf(myRank, "wait\n");
            ring_rx_chunk(nbArrEle, nbProc, myRank, buff, &tag);
            mpi_printf(myRank, "rxed1\n");

            if(tag == myRank){
                proc_my_chunk(nbArrEle, nbProc, buff);
                
                ring_tx(myRank, nbProc, (char*)buff, sizeof(buff), 0);
                mpi_printf(myRank, "protxed\n");
                continue;
            }

            ring_tx(myRank, nbProc, (char*)buff, sizeof(buff), tag);
            mpi_printf(myRank, "txed\n");
        }
    }
}

void broadcast_ring(int myRank, int nbArrEle, int nbProc){
    if(myRank != 0){
        int buff[nbArrEle/nbProc];
            int tag;
            mpi_printf(myRank, "wait\n");
            ring_rx_chunk(nbArrEle, nbProc, myRank, buff, &tag);
            mpi_printf(myRank, "rxed1\n");

            if(tag == myRank){
                proc_my_chunk(nbArrEle, nbProc, buff);
                
                ring_tx(myRank, nbProc, (char*)buff, sizeof(buff), 0);
                mpi_printf(myRank, "protxed\n");
            }
    } else {
    }
    // send 1st to everyone
    // then do my computation
    // wait for result from everyone
    int arr[nbArrEle]; // simple case of nbProc = array size
    for (size_t i = 0; i < NB_ELE(arr); i++){
        arr[i] = i;
    }
    mpi_printf(myRank, "inited\n");
    for(int gi = 1; gi < nbProc; ++gi){ // cause 0 is txer
        ring_tx(myRank, nbProc, (char*)&arr[gi * NB_ELE(arr)/nbProc], NB_ELE(arr)/nbProc * sizeof(int), gi);
    }
    mpi_printf(myRank, "txed\n");
    proc_my_chunk(nbArrEle, nbProc, arr);
    mpi_printf(myRank, "proced\n");

    for(int gi = 1; gi < nbProc; ++gi){ // cause 0 is rxer
        int buff[NB_ELE(arr)/nbProc];
        int tag;
        mpi_printf(myRank, "wait\n");
        ring_rx_chunk(nbArrEle, nbProc, myRank, buff, &tag);
        assert(tag == 0);
        mpi_printf(myRank, "rxed1\n");
        int arrIdx = tag * NB_ELE(arr)/nbProc;
        memcpy((char*)&arr[arrIdx], buff, sizeof(buff));
    }
    mpi_printf(myRank, "Array output:");
    for (size_t i = 0; i < NB_ELE(arr); i++){
        printf("%d ", arr[i]);
    }
    printf("\n");

}

int main(int argc, char **argv){

    MPI_Init(&argc, &argv);

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    broadcast_ring(rank, size, size);

    MPI_Finalize();

    return 0;
}