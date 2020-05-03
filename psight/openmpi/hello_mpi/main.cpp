//
// Created by csk on 03/05/2020.
//
#include<stdio.h>
#include<mpi.h>
#include<math.h>

#define ELE_SIZ(x) sizeof(x[0])
#define NB_ELE(x) sizeof(x)/ELE_SIZ(x)

using namespace std;
int main(int argc, char **argv){

    MPI_Init(&argc, &argv);

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size); // nb processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(rank == 0){
        cout << "Enter a number " << endl;
        int n;
        cin >> n;
        for (int i = 1; i < size; ++i) { // starts with 1 coz, this is being sent from 0
            // buff, nb ele, type(mpi macro), dest, tag(will see later), comm obj(proc handler)
            MPI_Send(&n, 1, MPI_INT, i, 0, MPI_COMM_WORLD); // sending data to everyone else
        }
    } else {
        int n;
        // buff, nb ele, type, source, tag, comm obj, status
        MPI_Recv(&n, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        cout << n << "^" << rank << "=" << pow(n, rank) << endl;
    }


    MPI_Finalize();

    return 0;
}