
#include  <stdio.h>
#include  <stdlib.h>
#include  <sys/types.h>
#include  <unistd.h>

#include <mpi.h>
#include <pthread.h>

#include "svipc.h"

int spawn (char* program, char** arg_list)
{
  pid_t child_pid;
  /* Duplique ce processus. */
  child_pid = fork ();
  if (child_pid != 0)
     /* Nous sommes dans le processus parent. */
     return child_pid;
  else {
     /* Exécute PROGRAM en le recherchant dans le path. */
     execvp (program, arg_list);
     /* On ne sort de la fonction execvp uniquement si une erreur survient. */
     fprintf (stderr, "une erreur est survenue au sein de execvp\n");
     abort ();
  }
}
void *tSystem(void *arg)
{
  char* arg_list=(char*)arg;
  system(arg_list);
  return (void*)NULL;
}
int main (int argc, char *argv[])
{
  //int provided;
  //MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &provided);
  //printf ("provided %d\n",provided);
  MPI_Init(&argc, &argv);

  int my_rank;
  int size;
  /* DETERMINE RANK OF THIS PROCESSOR*/
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  /*DETERMINE TOTAL NUMBER OF PROCESSORS*/
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if(my_rank==0) {
  int shm_key = 1234;
  int sem_key = 2345;
  char id[]="donnee";
  int i=0;

  int status = svipc_shm_init(shm_key, 10);
  status = svipc_sem_init(sem_key, 10);

  slot_array arr;
  arr.countdims = 2;
  arr.number = (int *)malloc(arr.countdims * sizeof(*arr.number));
  arr.number[0]=100;
  arr.number[1]=100;
  arr.typeid = SVIPC_LONG;

  arr.data=malloc(100*100*sizeof(long));
  long *d=(long*)arr.data;

  for(i=0; i<10000; i++){
    d[i]=i;
  }
  status = svipc_shm_write(shm_key, id, &arr, 0);

//  char* arg_list[] = {
//     "yorick",     /* argv[0], le nom du programme. */
//     "-batch",
//     "test.i",
//     NULL      /* La liste d'arguments doit se terminer par NULL.  */
//  };
//  spawn ("yorick", arg_list);

  char* arg_list[] = {
     "yorick",     /* argv[0], le nom du programme. */
     "-batch",
     "test.i",
     NULL      /* La liste d'arguments doit se terminer par NULL.  */
  };

  // intercommunicator
  MPI_Comm everyone;
  MPI_Comm_spawn("mpi_yorick", arg_list, 1,
           MPI_INFO_NULL, 0, MPI_COMM_SELF, &everyone,
           MPI_ERRCODES_IGNORE);

  printf ("on attend Yorick\n");
  status = svipc_semtake(sem_key, 0, 1, 20);
  status = svipc_shm_read(shm_key, id, &arr, 0);

  MPI_Send(d, 10000, MPI_LONG, 1, 1, MPI_COMM_WORLD);

  free(arr.number);
  free(arr.data);

  status = svipc_shm_cleanup(shm_key);
  status = svipc_sem_cleanup(sem_key);
  } else {
    MPI_Status stat;
    long *d=malloc(100*100*sizeof(long));
    MPI_Recv(d, 10000, MPI_LONG, 0, 1, MPI_COMM_WORLD, &stat);
    int i;
    for(i=0; i<10; i++){
      printf("%d ",d[i]);
    }
    printf("\n");
    free(d);
  }
  MPI_Finalize();
  return 0;
}
