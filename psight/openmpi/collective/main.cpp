#include<stdio.h>
#include<mpi.h>
#include<math.h>
#include<time.h>
#include<stdlib.h>

#define ELE_SIZ(x) sizeof(x[0])
#define NB_ELE(x) sizeof(x)/ELE_SIZ(x)

void mpi_printf(int rank, const char *fmt, ...){
    printf("%d : ", rank);
    va_list args;
    va_start(args, fmt);
    vprintf(fmt, args);
    va_end(args);
    fflush(stdout);
    return;
}

void mpi_bcast(char *buff, int cnt){
    // 0 -> txer
    MPI_Bcast(buff, cnt, MPI_CHAR, 0, MPI_COMM_WORLD);
}


int broadcast(int argc, char **argv){

    MPI_Init(&argc, &argv);

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size); // nb processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int value = 45;
    mpi_bcast((char*)&value, sizeof(value)); // same api with same signature to both tx and rx, pretty cool ha!
    
    if(rank == 0){
        mpi_printf(rank, "bcasted %d\n", value);   
    } else {
        mpi_printf(rank, "rxed %d from 0\n", value);    
    }
    
    MPI_Barrier(MPI_COMM_WORLD); // to understand, comment this n see

    mpi_printf(rank, "buhbye!\n");

    MPI_Finalize();

    return 0;
}

int reduce(int argc, char **argv){

    MPI_Init(&argc, &argv);

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size); // nb processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    const int itemsPerProcess = 10;
    const int count = size * itemsPerProcess;

    int *data = (int*)malloc(count * sizeof(count));

    if(rank == 0){
        mpi_printf(rank, "numbers generated are: ");
        for (size_t i = 0; i < count; i++){
            data[i] = rand() % 10;
            printf("%d ", data[i]);
        }
        printf("\n");
    }

    int *localData = (int*)malloc(itemsPerProcess * sizeof(int));

    // array to be scattered, nb items/proc to be scattered, data type, where to scatter the data to, how many to be recieved, data type to be recieved, scatterer, proc obj
    MPI_Scatter(data, itemsPerProcess, MPI_INT, localData, itemsPerProcess, MPI_INT, 0, MPI_COMM_WORLD); // 0 -> rank of whoever's 'scattering'
    free(data);
        
    int localSum = 0;
    for (size_t i = 0; i < itemsPerProcess; i++){
        localSum += localData[i];
    }
    free(localData);

    mpi_printf(rank, "localsum: %d\n", localSum);
    
    int globalSum;
    // local var to be reduced, the var to which reduced output is written,  data type, operation to be performed, scatterer, com/proc obj
    MPI_Reduce(&localSum, &globalSum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if(rank == 0){
        mpi_printf(rank, "global sum: %d\n", globalSum);
    }

    MPI_Finalize();

    return 0;
}

int main(int argc, char **argv){
    // return broadcast(argc, argv);
    srand(time(NULL));
    return reduce(argc, argv);
}