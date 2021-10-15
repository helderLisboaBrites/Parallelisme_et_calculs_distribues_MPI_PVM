/*
*    tokenring example using PVM 3 and group functions
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include "pvm3.h"

main(int argc, char ** argv)
{
    int NPROC = 8;		/* default nb of proc */
    int mytid;                  /* my task id */
    int *tids;                  /* array of task id */
    int me;                     /* my process number */
    int i;

    /* enroll in pvm */
    mytid = pvm_mytid();

    /* Join a group and if I am the only instance     */
    /* i.e. pvm_gsize = 0 spawn more copies of myself */

    /* gathering of tids in the tids array is optional */
    tids = (int *)malloc ( sizeof(int) * NPROC);

    me = pvm_joingroup( "foo" ); 
    printf("P%d (%x)\n", me, mytid); fflush (stdout);

    if( pvm_gsize("foo") == 1 ) { 
      tids[0] = pvm_mytid();
      for( i = 1; i < NPROC; i++) {
	pvm_spawn( argv[0], (char**)0, 0, "", 1, &tids[i]);
	printf("P%d : spawning task %d tid : %x\n", me, i, tids[i]); fflush (stdout);
      }
      /* Wait for everyone to startup before proceeding. */
      pvm_barrier( "foo", NPROC );
    }
    else { /* We're several : I am a spawned process */
      /* Wait for everyone to startup before proceeding. */
      pvm_barrier( "foo", NPROC );

      NPROC = pvm_gsize("foo");
      printf("P%d : nproc: %d\n", me, NPROC); fflush (stdout);
    }
    pvm_freezegroup ( "foo", NPROC );
    for ( i = 0; i < NPROC; i++) tids[i] = pvm_gettid ( "foo", i);

/*--------------------------------------------------------------------------*/

    /* all the tasks are equivalent at that point */

    { char sortie [255];
    sprintf ( sortie, "P%d (%x) tids :", me, mytid );
    for ( i = 0; i < NPROC; sprintf ( sortie+strlen(sortie), " %x", tids[i++]) );
    printf ( "%s\n", sortie);
    fflush (stdout);
    }
     
     dowork( me, NPROC );

     /* program finished leave group and exit pvm */
     pvm_lvgroup( "foo" );
     printf ("P%d (%x) exiting\n", me, mytid); fflush ( stdout );
     pvm_exit();
     exit(1);
}

/* Simple example passes a token around a ring */

dowork( me, nproc )
     int me;
     int nproc;
{
     int token;
     int src, dest;
     int count  = 1;
     int stride = 1;
     int msgtag = 4;

     /* Determine neighbors in the ring */
     src = pvm_gettid ("foo", (me-1+nproc) % nproc );
     dest= pvm_gettid ("foo", (me+1)%nproc );

     if( me == 0 )
     { 
        token = dest;
        pvm_initsend( PvmDataDefault );
        pvm_pkint( &token, count, stride );
        pvm_send( dest, msgtag );
        pvm_recv( src, msgtag );
        printf("token ring done\n");
     }
     else
     {
        pvm_recv( src, msgtag );
        pvm_upkint( &token, count, stride );
        pvm_initsend( PvmDataDefault );
        pvm_pkint( &token, count, stride );
        pvm_send( dest, msgtag );
     }
}
