#include "mpi.h"
#include <stdio.h>
#include <sys/types.h>

#include <stdlib.h>

#include <sys/time.h>
#include <time.h>
#include <string.h>
#include <errno.h>


double f(double x){
    return 1.0f/(1+x*x);
}

double integrale(int N, int me, int nproc){
    double aire;
    for(int k = me*N/nproc; k < (me+1)*N/nproc; k++){
        aire += (1.0f/N) * f(k*1.0f/N);
    }
    return aire;
}



int main(int argc, char *argv[]){

    int me;
    MPI_Status status;
    int N;
    int nproc;
    double trueResultat;
    if (argc != 2)
    {
        printf("Usage: %s <N> \n", argv[0]);
        exit(-1);
    }
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    MPI_Comm_size( MPI_COMM_WORLD, &nproc );

    N = atoi(argv[1]);
    MPI_Bcast(&N,1,MPI_INT,0,MPI_COMM_WORLD);


    double resultat = integrale(N,me,nproc);
    printf("resultat aire = %lf moi = %d N=%d\n", resultat, me,N);
    MPI_Reduce(&resultat, &trueResultat, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD);

    if(me == 0){
        printf("Resultat de pi = %lf", trueResultat*4);
    }
    MPI_Finalize();

    return 0;
}