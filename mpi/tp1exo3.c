#include "mpi.h"
#include <stdio.h>
#include <sys/types.h>

#include <stdlib.h> 
#include <sys/time.h>
#include <time.h>
#include <string.h>
#include <errno.h>

double f(double x)
{
    return 1.0f / (1 + x * x);
}

double integrale(int N, int me, int nproc)
{
    double aire;
    for (int k = me * N / nproc; k < (me + 1) * N / nproc; k++)
    {
        aire += (1.0f / N) * f(k * 1.0f / N);
    }
    return aire;
}

void sendLigne(int me, int nproc, int N, double *tab)
{
    if (me == 0)
    {
        MPI_Isend(tab + N * (nproc-1), N, MPI_DOUBLE_PRECISION, 1, -1, MPI_COMM_WORLD);
    }
    else if (me == nproc-1)
    {
        MPI_Isend(tab, N, MPI_DOUBLE_PRECISION, nproc - 1, -1, MPI_COMM_WORLD);
    }
    else
    {
        MPI_Isend(tab, N, MPI_DOUBLE_PRECISION, me - 1, -1, MPI_COMM_WORLD);
        MPI_Isend(tab + N * (nproc-1), N, MPI_DOUBLE_PRECISION, me + 1, -1, MPI_COMM_WORLD);
    }
}

void recvLigne(int me, int nproc, int N, double *tab, double *ligneDown, double *ligneUp,MPI_Status *status)
{
    if (me == 0)

    {
        MPI_Recv(ligneDown, N, MPI_DOUBLE_PRECISION, 1, -1, MPI_COMM_WORLD,status);
    }
    else if (me == nproc-1)
    {
        MPI_Recv(ligneUp, N, MPI_DOUBLE_PRECISION, nproc - 1, -1, MPI_COMM_WORLD,status);
    }
    else
    {
        MPI_Recv(ligneUp, N, MPI_DOUBLE_PRECISION, me - 1, -1, MPI_COMM_WORLD,status);
        MPI_Recv(ligneDown, N, MPI_DOUBLE_PRECISION, me + 1, -1, MPI_COMM_WORLD,status);
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

    double *ligneUp = (double *)calloc(N, sizeof(double));
    double *ligneDown = (double *)calloc(N, sizeof(double));

    double *tab = (double *) calloc(N * N / nproc , sizeof(double));

    printf("azeazeazeazeaze\n"); fflush(stdout);
    for (int i = 0; i < nproc; i++)
    {
        for (int j = 0; j < N; j++)
        {
            *(tab + i * N + j) = me;
        }
    }

    for (int i = 0; i < N; i++)
        *(ligneUp + i) = -1;
    for (int i = 0; i < N; i++)
        *(ligneDown + i) = -1;

    MPI_Barrier(MPI_COMM_WORLD);
    sendLigne(me, nproc, N, tab);

printf("11111111111111111111111\n"); fflush(stdout);
    recvLigne(me, nproc, N, tab, ligneDown, ligneUp,&status);

    for (int n = 0; n < nproc; n++)
    {
        if (me == n)
        {
            for (int u = 0; u < N; u++)
            {
                printf("%lf", ligneUp + u);
            }
            printf("\n");
            for (int i = 0; i < nproc; i++)
            {

                for (int j = 0; j < N; j++)
                {
                    printf("%lf", tab + i * j);
                }
                printf("\n");
            }
            for (int v = 0; v < N; v++)
            {
                printf("%lf", ligneUp + v);
            }
            printf("\n");
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    return 0;
}