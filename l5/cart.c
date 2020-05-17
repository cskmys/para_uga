#include "cart.h"

#include <string.h>
#include <stdlib.h>

#include "utils.h"

void createCartTopology(int row, int col, MPI_Comm *comm2d){
    int nbNodes = row * col;
    int dim[NB_DIM];
    dim[0] = row;
    dim[1] = col;
    MPI_Dims_create(nbNodes, NB_DIM, dim);
    int wrapArnd[NB_DIM];
    memset(wrapArnd, 0, sizeof(wrapArnd));
    int reorder = 1;
    int err;
    if((err = MPI_Cart_create(MPI_COMM_WORLD, NB_DIM, dim, wrapArnd, reorder, comm2d)) != 0){
        mpi_printf("err: %d creating cart\n", err);
        exit(-1);
    }
    return;
}

void getCartCoord(MPI_Comm comm2d,  int *coord2d){
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Cart_coords(comm2d, rank, NB_DIM, coord2d);
    return;
}

int getCartRank(MPI_Comm comm2d){
    int coord2d[NB_DIM];
    getCartCoord(comm2d, coord2d);
    int cartRank;
    MPI_Cart_rank(comm2d, coord2d, &cartRank);
    return cartRank;
}

void splitCartByRow(MPI_Comm comm2d, MPI_Comm *comm1dRow){
    int freeCoords[NB_DIM] = { 0, 1};
    MPI_Cart_sub(comm2d, freeCoords, comm1dRow);
}

void getRowCoord(MPI_Comm comm1dRow,  int *coord1d){
    int rowRank = getRowRank(comm1dRow);
    MPI_Cart_coords(comm1dRow, rowRank, 1, coord1d);
    return;
}

int getRowRank(MPI_Comm comm1dRow){
    int rowRank;
    MPI_Comm_rank(comm1dRow, &rowRank);
    return rowRank;
}

void splitCartByCol(MPI_Comm comm2d, MPI_Comm *comm1dCol){
    int freeCoords[NB_DIM] = { 1, 0};
    MPI_Cart_sub(comm2d, freeCoords, comm1dCol);
}

void getColCoord(MPI_Comm comm1dCol,  int *coord1d){
    int colRank = getColRank(comm1dCol);
    MPI_Cart_coords(comm1dCol, colRank, 1, coord1d);
    return;
}

int getColRank(MPI_Comm comm1dCol){
    int colRank;
    MPI_Comm_rank(comm1dCol, &colRank);
    return colRank;
}

void delCol(MPI_Comm *comm1dCol){
    MPI_Comm_free(comm1dCol);
}

void delRow(MPI_Comm *comm1dRow){
    MPI_Comm_free(comm1dRow);
}

void delCartTopology(MPI_Comm *comm2d){
    MPI_Comm_free(comm2d);
}