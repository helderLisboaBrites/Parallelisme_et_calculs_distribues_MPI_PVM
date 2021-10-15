
/*  tokenring example using PVM 3.4 
    - uses sibling() to determine the nb of spawned tasks (xpvm and pvm> ok)
    - uses group for token ring communication
*/

#include <stdio.h>
#include <sys/types.h>
#include "pvm3.h"

#include <stdlib.h>

#include <sys/time.h>
#include <time.h>
#include <string.h>
#include <errno.h>

#define GRPNAME "gauss"

void matrix_display ( double *tab,int  N ) {
  int i,j;

  for (i=0; i<N; i++) {
    for (j=0; j<N; j++) {
      printf ("%8.2f ", *(tab+i*N+j) );fflush(stdout);
    }
    printf ("\n");fflush(stdout);
  }

}

void gauss ( double * tab, int N,int nproc,int me ) {
  int i,j,k;
  double pivot;

  for ( k=0; k<N-1; k++ ){ /* mise a 0 de la col. k */
    /* printf (". "); */
    if ( fabs(*(tab+k+k*N)) <= 1.0e-11 ) {
      printf ("ATTENTION: pivot %d presque nul: %g\n", k, *(tab+k+k*N) );
      exit (-1);
    }
    for ( i=(int)(k+1)/me; i<N/nproc; i++ ){ /* update lines (k+1) to (n-1) */
      pivot = - *(tab+k+i*N/nproc) / *(tab+k+k*N/nproc);
      for ( j=k; j<N; j++ ){ /* update elts (k) - (N-1) of line i */
	*(tab+j+i*N) = *(tab+j+i*N) + pivot * *(tab+j+k*N);
      }
      /* *(tab+k+i*N) = 0.0; */
    }
  }
  printf ("\n");
}


void matrix_load(char nom[], double *tab, int N, int me, int nproc)
{
  FILE *f;
  int i, j;
  double ligne[N];
  if (me == 0)
  {
    if ((f = fopen(nom, "r")) == NULL)
    {
      perror("matrix_load : fopen ");
    }

    for (i = 0; i < N; i++)
    {
      for (j = 0; j < N; j++)
      {
        fscanf(f, "%lf", ligne + j);
        //printf("%lf \n", ligne[j]);fflush(stdout);
      }
      int ilocal = (int)i / nproc;
      if (i % nproc == 0)
      {
        memcpy(tab + ilocal * N, ligne, N * sizeof(double));
      }
      else
      {

        int dest = pvm_gettid(GRPNAME, i % nproc);
        printf("send %d \n", dest);fflush(stdout);
        pvm_initsend( PvmDataDefault );
        pvm_pkdouble(ligne, N, 1);
        pvm_send(dest, me);
      }
    }
    fclose(f);
  }
  else
  {
    for (int i = 0; i < N / nproc; i++)
    {
      pvm_recv(-1, -1);
      pvm_upkdouble(tab + i*N, N, 1);
    }
  }
}
/*spawn -> -4 tpgauss 16 m16*/
void matrix_save(char nom[], double *tab, int N, int me, int nproc)
{
  FILE *f;
  int i, j;
  double ligne[N];
printf("azeazeazeaze %d \n", me);fflush(stdout);
  if (me == 0)
  {
    if ((f = fopen(nom, "w")) == NULL)
    {
      perror("matrix_save : fopen ");
    }
    
    //On recupere les lignes du tab
    for (int i = 0; i < N ; i++)
    {
      if(i%nproc == 0){
        if(i==0){
          printf("bonjour\n");fflush(stdout);
          matrix_display(tab,N);
        }
        for (int j = 0; j < N ; j++)
          fprintf(f, "%8.2f ", *(tab+i*N + j));
        fprintf(f, "\n");
      }else{
        printf("reciv %d \n", me);fflush(stdout);
        pvm_recv(-1, i% nproc);
        pvm_upkdouble(ligne, N, 1);
        for (int j = 0; j < N ; j++)
          fprintf(f, "%8.2f ", ligne[j]);
        fprintf(f, "\n");
        
      }
    }
    fclose(f);
  }
  else
  {
    for (int i = 0; i < N / nproc ; i++)
    {
      int dest = pvm_gettid(GRPNAME, 0);
      pvm_initsend( PvmDataDefault );
      pvm_pkdouble(tab+i*N, N, 1);
      pvm_send(dest, me);
    }
  }
}

main(int argc, char **argv)
{
  int NPROC = 8; /* default nb of proc */
  int mytid;     /* my task id */
  int *tids;     /* array of task id */
  int me;        /* my process number */
  int i;
  int N;
  char nom[255];

  if (argc != 3)
  {
    printf("Usage: %s <matrix size> <matrix name>\n", argv[0]);
    exit(-1);
  }
  N = atoi(argv[1]);
  strcpy(nom, argv[2]);

  /* enroll in pvm */
  mytid = pvm_mytid();

  /* determine the size of my sibling list */
  NPROC = pvm_siblings(&tids);
  /* WARNING: tids are in order of spawning, which is different from
       the task index JOINING the group */

  me = pvm_joingroup(GRPNAME); /* me: task index in the group */
  pvm_barrier(GRPNAME, NPROC);
  pvm_freezegroup(GRPNAME, NPROC);
  for (i = 0; i < NPROC; i++)
    tids[i] = pvm_gettid(GRPNAME, i);

  /*--------------------------------------------------------------------------*/
  /*           all the tasks are equivalent at that point                     */

  dowork(me, tids, NPROC, nom, N);

  pvm_lvgroup(GRPNAME);
  pvm_exit();
}

/* Simple example passes a token around a ring */

dowork(int me, int tids[], int nproc, char nom[], int N)
{
  int i, j, k;
  double *tab, pivot;
  FILE *f;
  struct timeval tv1, tv2; /* for timing */
  int duree;

  /* déclaration du tableau de taille N*N/nproc pour chaque processeur
  */
  if ((tab = malloc(N * N / nproc * sizeof(double))) == NULL)
  {
    printf("Cant malloc %d bytes\n", (int)(N * N / nproc) * sizeof(double));
    exit(-1);
  }

  gettimeofday(&tv1, (struct timezone *)0);
  matrix_load(nom, tab, N, me, nproc);

  gettimeofday(&tv2, (struct timezone *)0);
  duree = (tv2.tv_sec - tv1.tv_sec) * 1000000 + tv2.tv_usec - tv1.tv_usec;
  printf("loading time: %10.8f sec.\n", duree / 1000000.0);
  fflush(stdout);
  //matrix_display(tab,N);
  
  gettimeofday(&tv1, (struct timezone *)0);
  gauss(tab, N,nproc,me);
  gettimeofday(&tv2, (struct timezone *)0);
  duree = (tv2.tv_sec - tv1.tv_sec) * 1000000 + tv2.tv_usec - tv1.tv_usec;
  printf("computation time: %10.8f sec.\n", duree / 1000000.0);

  sprintf(nom + strlen(nom), ".result");
  matrix_save(nom, tab, N, me, nproc);
}
