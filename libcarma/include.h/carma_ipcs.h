#ifndef CARMA_IPCS_H
#define CARMA_IPCS_H
//cuda
#include <cuda.h>
//shms
#include <sys/mman.h>
#include <sys/stat.h>        /* For mode constants */
#include <fcntl.h>           /* For O_* constants */
#include <semaphore.h>
//malloc
#include <stdlib.h>
//getpid
#include <sys/types.h>
#include <unistd.h>

#include <cstring>
#include <stdio.h>
#include <errno.h>
#include <limits.h>

#include <map>

#define STAMP(fmt, args...) fprintf(stderr, "[%s@%d]:" fmt, __FUNCTION__, __LINE__, ## args)

//todo: CIPCS_CHECK macro




//sizeof(CUipcMemHandle) == sizeof(CUipcEventHandle)
//interessant !!
//possibilité d'unifier le tout et de déclarer moins de struct
//a creuseer


typedef struct sh_cu_dptr_st{
  char              name[NAME_MAX+1];
  unsigned int      nb_proc;
  pid_t             owner;
  CUipcMemHandle    handle;
  sem_t             var_mutex;
}sh_dptr;

typedef struct sh_cu_event_st{
  char              name[NAME_MAX+1];
  unsigned int      nb_proc;
  pid_t             owner;
  CUipcEventHandle  handle;
  sem_t             var_mutex;
}sh_event;

typedef struct sh_buffer_st{
  char              name[NAME_MAX+1];
  size_t            size;
  unsigned int      nb_proc;
  void              *p_shm;
  sem_t             mutex;
}sh_buffer;

typedef struct sh_barrier_st{
  char              name[NAME_MAX+1];
  bool              valid;
  unsigned int      proc_cnt;
  unsigned int      val;
  sem_t             b_sem;
  sem_t             var_mutex;
}sh_barrier;


class carma_ipcs{
private:
  std::map<unsigned int, sh_dptr *> dptrs;
  std::map<unsigned int, sh_event *> events;
  std::map<unsigned int, sh_barrier *> barriers;
  std::map<unsigned int, sh_buffer> buffers;

  void * create_shm(const char *name, size_t size);
  void * get_shm(const char *name);
  void free_shm(const char *name, void *p_shm, size_t size);
  void close_shm(void *p_shm, size_t size);
  void complete_clean();

public:
  /*
   general purpose methods
  */
  carma_ipcs();
  ~carma_ipcs();


  /*
    Cuda handles methods
  */
  //register a device pointer, id must be a non nul argument
  int register_cudptr(CUdeviceptr dptr, unsigned int id);
  //register an event, id must be a non nul argument
  int register_cuevent(CUevent event, unsigned int id);

  //get a memory handle
  int get_memHandle(unsigned int id, CUipcMemHandle *phandle);
  //get a event handle
  int get_eventHandle(unsigned int id, CUipcEventHandle *phandle);

  //free a memory handle shared mem space
  void free_memHandle(unsigned int id);
  //free a event handle shared event space
  void free_eventHandle(unsigned int id);


  /*
    Transfer via CPU memory methods
  */
  //allocation of the shm for memory tranfers
  int alloc_memtransfer_shm(size_t bsize, unsigned int id, void *shm);
  //get tranfer shm
  int get_memtransfer_shm(unsigned int id, void *shm);
  //free transfer shm ref by id
  void free_memtransfer_shms(unsigned int id);


  /*
    Barrier methods
  */
  int init_barrier(unsigned int id, unsigned int value);

  int wait_barrier(unsigned int id);

  void free_barrier(unsigned int id);

};

#endif //CARMA_IPCS_H
