//
// Created by csk on 03/05/2020.
//
#include<stdio.h>
#include<string.h>
#include<mpi.h>
#include<math.h>
#include <stdarg.h>

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

typedef struct{
    char name[16];
    int age;

    void print(int rank){
        mpi_printf(rank, "%s is %d years old", name, age);
    }
    void setName(char *str){
        strncpy(name, "csk", strnlen(str, sizeof(name)) + 1);
    }
    void setAge(int val){
        age = val;
    }
}Person;

typedef union{ // not always recommended due to potential endianess difference over machines across network
    Person person;
    char buff[sizeof(Person)];
}UPerson;

int main(int argc, char **argv){

    MPI_Init(&argc, &argv);

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size); // nb processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(rank == 0){
        UPerson p;
        p.person.setName((char*)"csk");
        p.person.setAge(27);
        for (int i = 1; i < size; ++i) {
            MPI_Send(p.buff, sizeof(p.buff), MPI_CHAR, i, 0, MPI_COMM_WORLD); // non-blocking, doesnt care if someone received or not
        }
    } else {
        // we dont know in advance how much to receive
        MPI_Status status;
        MPI_Probe(0, 0, MPI_COMM_WORLD, &status); // blocking wait for msg from rank 0
        int count;
        MPI_Get_count(&status, MPI_CHAR, &count); // using char as type for size of msg i,e, max 255 bytes of msg at once??
        mpi_printf(rank, "have %d bytes\n", count);
        UPerson p;
        MPI_Recv(p.buff, count, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status); // blocking call
        mpi_printf(rank, "got %d bytes\n", count);
        p.person.print(rank);
        printf("\n");
    }

    MPI_Finalize();

    return 0;
}