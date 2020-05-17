#include<stdio.h>
#include<stdarg.h>
#include<string.h>
#include<stdlib.h>

#include<mpi.h>

#include "utils.h"

#ifdef DBG_ON_FD

FILE *pFd[MAX_NB_PROC] = {NULL};

void init_comm_pipe(int rank){
    char *filNam = strdup("f000.log");
    snprintf(filNam, strlen(filNam), "f%02d.log", rank);

    if((pFd[rank] = fopen(filNam, "wb")) < 0){
        printf("%d: open failed!\n", rank);
        free(filNam);
        while(1);
    }
    free(filNam);
}

void deinit_comm_pipe(int rank){
    fclose(pFd[rank]);
}
#endif

void mpi_fprintf(FILE *f, const char *fmt, va_list args){
    
    vfprintf(f, fmt, args);
    
    fflush(f);
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

void mpi_printf(const char *fmt, ...){
#ifdef DBG_ON_FD
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    va_list args;
    va_start(args, fmt);

    mpi_fprintf(pFd[rank], fmt, args);

    va_end(args);
#endif
    return;
}

int getIdx(int nbEle, int idx){
    idx = idx % nbEle;
    if(idx < 0){
        idx = nbEle + idx;
    }
    return idx;
}

double experiments [NBEXPERIMENTS] ;
double average_time (void){
    unsigned int i ;
    double s = 0 ;

    for (i = 0; i < NBEXPERIMENTS; i++){
        s = s + experiments [i] ;
    }
    return s / NBEXPERIMENTS;
}

void set_exp_time(int trial, double tim){
    trial = trial % NB_ELE(experiments);
    experiments[trial] = tim;
}