/*
 * mpi_yorick.c
 *
 *  Created on: Apr 9, 2014
 *      Author: sevin
 */

#include  <stdio.h>
#include  <stdlib.h>
#include  <sys/types.h>
#include  <unistd.h>
#include <mpi.h>

//#include "svipc.h"

int main (int argc, char *argv[])
{
  MPI_Init(&argc, &argv);
  if(argc==4){
    char cmd[128];
    sprintf(cmd, "%s %s %s", argv[1], argv[2], argv[3]);
    system(cmd);//, &argv[1]);
    printf("j'ai quitté Yorick\n");
  } else {
    fprintf(stderr,"Parametres incorrects\n");
  }
  MPI_Finalize();
  return EXIT_SUCCESS;
}


