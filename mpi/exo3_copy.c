#include "mpi.h"
#include <stdio.h>
#include <sys/types.h>

#include <stdlib.h> 
#include <sys/time.h>
#include <time.h>
#include <string.h>
#include <errno.h>


void sendLigne(int me, int nproc, int N, int totalLignes, double *tab)
{
    MPI_Request r;
    if (me == 0)
    {
        MPI_Isend(tab + N * (totalLignes-2), N, MPI_DOUBLE_PRECISION, 1, 0, MPI_COMM_WORLD,&r);
    }
    else if (me == nproc-1)
    {
        MPI_Isend(tab+N, N, MPI_DOUBLE_PRECISION, nproc - 2, 0, MPI_COMM_WORLD,&r);
    }
    else
    {
        MPI_Isend(tab+N, N, MPI_DOUBLE_PRECISION, me - 1, 0, MPI_COMM_WORLD,&r);
        MPI_Isend(tab + N * (totalLignes-2), N, MPI_DOUBLE_PRECISION, me + 1, 0, MPI_COMM_WORLD,&r);
    }
}

void recvLigne(int me, int nproc, int N, double *tab, int totalLignes,MPI_Status *status)
{
    if (me == 0)

    {
        MPI_Recv(tab + N * (totalLignes-1), N, MPI_DOUBLE_PRECISION, 1, MPI_ANY_TAG, MPI_COMM_WORLD,status);
    }
    else if (me == nproc-1)
    {
        MPI_Recv(tab, N, MPI_DOUBLE_PRECISION, nproc - 2, MPI_ANY_TAG, MPI_COMM_WORLD,status);
    }
    else
    {
        MPI_Recv(tab, N, MPI_DOUBLE_PRECISION, me - 1, MPI_ANY_TAG, MPI_COMM_WORLD,status);
        MPI_Recv(tab + N * (totalLignes-1), N, MPI_DOUBLE_PRECISION, me + 1, MPI_ANY_TAG, MPI_COMM_WORLD,status);
    }
}


int main(int argc, char *argv[])
{

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
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    N = atoi(argv[1]);

    int totalLignes = N/nproc +2; 

    double *tab = (double *) calloc(N * totalLignes , sizeof(double)); 
    // plus deux ligne celle d'avant et celle d'apres 
    // donc la ligne 0 est la ligne d'avant et la ligne N/nproc+1 la ligne d'apres

    //printf("totalLignes %d\n",totalLignes);fflush(stdout);
    for (int i = 0; i < totalLignes; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if(i == 0 || i == totalLignes-1 )
                *(tab + i * N + j) = -1;
            else 
                *(tab + i * N + j) = me;
        }
    }
    /*
    for(int k =0; k < nproc; k++){
        if(k == me)
            for(int i = 0; i<totalLignes; i++){
                for(int j = 0; j<N; j++){
                    printf(" %lf ",*(tab + i * N + j));fflush(stdout);
                }
                printf("\n");fflush(stdout);
                }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    */
    
    MPI_Barrier(MPI_COMM_WORLD);
    sendLigne(me, nproc, N, totalLignes,tab);

    recvLigne(me, nproc, N, tab, totalLignes, &status);

    MPI_Barrier(MPI_COMM_WORLD);
     for(int k =0; k < nproc; k++){
        if(k == me){
            for(int i = 0; i<totalLignes; i++){
                for(int j = 0; j<N; j++){
                    printf(" %lf ",*(tab + i * N + j));fflush(stdout);
                }
                printf("\n");fflush(stdout);
                }
        printf("\n");fflush(stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    
    free(tab);
   // free(ligneUp);
    return 0;
}