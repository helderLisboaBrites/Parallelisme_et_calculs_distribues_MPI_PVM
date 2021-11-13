#include "mpi.h"
#include <stdio.h>
#include <sys/types.h>

#include <stdlib.h> 
#include <sys/time.h>
#include <time.h>
#include <string.h>
#include <errno.h>
#include <math.h>

void sendLigne(int me, int nproc, int N, int totalLignes, double *tab)
{
    MPI_Request r;
    if (me == 0)
    {   //Le premier coin envoie à droite et en bas
        MPI_Isend(tab + N * (totalLignes-2), N, MPI_DOUBLE_PRECISION, 1, 0, MPI_COMM_WORLD,&r);
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


double convergence(double* f,double* fnew, int N, int totalLignes, int nproc){
    double err=0;
    double errResult;

    for(int i =1; i< totalLignes-1; i++){
        for(int j=0; j< N; j++){
            err += (*(fnew+i*N+j) - *(f + i*N + j) ) * (*(fnew+i*N+j) - *(f + i*N + j));
        }
    }
    MPI_Reduce(&err, &errResult, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(&errResult,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD);
    //printf("err : %lf \n",errResult);fflush(stdout);
    return sqrt( errResult );
}

void distributionLignes(int me,int nproc,int N, double* f, int totalLignes, MPI_Status *status){
    MPI_Barrier(MPI_COMM_WORLD);
    sendLigne(me, nproc, N, totalLignes,f);
    recvLigne(me, nproc, N, f, totalLignes, status);
    MPI_Barrier(MPI_COMM_WORLD);
}

void laplace(int me, double* f,double* fnew, int N, int totalLignes, int nproc, MPI_Status *status){
    int test =1;
    while(test){
        //printf("----------------------\n");fflush(stdout);

        for(int i =1; i< totalLignes-1; i++){
            for(int j = 1; j < N-1; j++){
                *(fnew + i*N + j) = 0.25 * ( *(f + (i+1)*N + j) + *(f + (i-1)*N + j)
                                           + *(f + i*N + j+1) + *(f + i*N + j-1));
                                           //printf(" %lf ",*(fnew + i*N + j));fflush(stdout);
            }
            //printf("\n");fflush(stdout);

        }
        test = (convergence(f, fnew, N, totalLignes, nproc) > 0.000001);
        
        for(int i =1; i< totalLignes-1; i++){
            for(int j = 0; j < N; j++){
                *(f+i*N+j) = *(fnew+i*N+j); 
            }
        }
        distributionLignes(me, nproc, N, f, totalLignes, status);
        
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



    int N_local = N/sqrt(nproc) +2;

    double *f = (double *) calloc(N_local * N_local , sizeof(double)); 
    double *fnew = (double *) calloc(N_local * N_local , sizeof(double)); 
    // plus deux ligne celle d'avant et celle d'apres 
    // donc la ligne 0 est la ligne d'avant et la ligne N/nproc+1 la ligne d'apres

    printf("totalLignes %d\n",N_local);fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
    for (int i = 0; i < N_local; i++)
    {
        for (int j = 0; j < N_local; j++)
        {
            if(i == 0 || i == N_local-1 || j==0 || j==N_local-1)
                *(f + i * N_local + j) = -1;
            else 
                *(f + i * N_local + j) = me;
        }
    }
    /*
    if(me ==3){
        for (int i = 0; i < N_local; i++)
        {
            for (int j = 0; j < N_local; j++){
                printf(" %2.3f ",*(f + i * N_local + j));fflush(stdout);
            }
        printf("\n");fflush(stdout);
        }
    printf("\n");fflush(stdout);
    }*/
/*
    //Partage des lignes up et down
    distributionLignes(me, nproc, N, f, totalLignes, &status);

    //Calcul de laplce
    laplace(me, f,fnew, N, totalLignes, nproc, &status);

*/
    MPI_Barrier(MPI_COMM_WORLD);
    //int i,fin;
/*
     for(int k =0; k < nproc; k++){
        if(k == me){
            if(me ==0) {i = 0; fin = totalLignes-1;}
            else if(me ==nproc-1) {i = 1; fin = totalLignes;}
            else {i = 1; fin = totalLignes-1;}
            
            for(i ; i<fin ; i++){
                for(int j = 0; j<N; j++){
                    printf(" %lf ",*(f + i * N + j));fflush(stdout);
                }
                printf("\n");fflush(stdout);
                }
        //printf("\n");fflush(stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
*/
/*
    for(int i=0; i < N_local*(int)sqrt(nproc) ;i++){
        for(int k =0; k < nproc; k++){
            if(k == me && (i/(int)sqrt(nproc)) ==me){
                for(int j=0;j<N_local;j++){
                    printf(" %.3f ",*(f + (i%N_local) * N_local + j));fflush(stdout);
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        if(me==0)
            printf("\n");fflush(stdout);
        MPI_Barrier(MPI_COMM_WORLD);
    }
*/
    MPI_Comm comm;
    MPI_Comm_split(MPI_COMM_WORLD, me / (int)sqrt(nproc), me, &comm);
    int row_rank, row_size;
    MPI_Comm_rank(comm, &row_rank);
    MPI_Comm_size(comm, &row_size);
    //printf(" moi = %d -- %d -- %d",me,row_rank,row_size);fflush(stdout);
    for(int l = 0; l < (int)sqrt(nproc); l++){
        //
        if( l * (int)sqrt(nproc) <= me && me < (l+1) * (int)sqrt(nproc) ){
            for(int i =0; i<N_local ; i++){
                for(int k = l * (int)sqrt(nproc); k <(l+1) * (int)sqrt(nproc); k++){
                    if(k==me){
                        for(int j =0; j<N_local ; j++){
                            printf(" %2.3f ",*(f + i * N_local + j));fflush(stdout);
                        }
                    }
                    MPI_Barrier(comm);
                }
                if(me == (l+1) * (int)sqrt(nproc)-1)
                    printf(" \n ");fflush(stdout);
                MPI_Barrier(comm);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    free(f);
    free(fnew);
    return 0;
}