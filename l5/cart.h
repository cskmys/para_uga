#ifndef CART_H
#define CART_H

#include <mpi.h>

#define NB_DIM 2

void createCartTopology(int row, int col, MPI_Comm *comm2d);
void getCartCoord(MPI_Comm comm2d,  int *coord2d);
int getCartRank(MPI_Comm comm2d);
void splitCartByRow(MPI_Comm comm2d, MPI_Comm *comm1dRow);
void getRowCoord(MPI_Comm comm1dRow,  int *coord1d);
int getRowRank(MPI_Comm comm1dRow);
void splitCartByCol(MPI_Comm comm2d, MPI_Comm *comm1dCol);
void getColCoord(MPI_Comm comm1dCol,  int *coord1d);
int getColRank(MPI_Comm comm1dCol);
void delCol(MPI_Comm *comm1dCol);
void delRow(MPI_Comm *comm1dRow);
void delCartTopology(MPI_Comm *comm2d);

#endif