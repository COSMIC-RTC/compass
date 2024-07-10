// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.
// -----------------------------------------------------------------------------

//! \file      carma_ipc.hpp
//! \ingroup   libcarma
//! \class     carma_ipc
//! \brief     this class provides the ipc features to CarmaObj
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#ifndef CARMA_IPCS_H
#define CARMA_IPCS_H
// cuda
#include <cuda.h>
// shms
#include <fcntl.h> /* For O_* constants */
#include <semaphore.h>
#include <sys/mman.h>
#include <sys/stat.h> /* For mode constants */
// malloc
#include <stdlib.h>
// getpid
#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>

#include <iterator>
#include <map>

#define STAMP(fmt, args...) \
  fprintf(stderr, "[%s@%d]:" fmt, __FUNCTION__, __LINE__, ##args)

// todo: CIPCS_CHECK macro

// sizeof(CUipcMemHandle) == sizeof(CUipcEventHandle)
// interessant !!
// possibilité d'unifier le tout et de déclarer moins de struct
// a creuseer

typedef struct sh_cu_dptr_st {
  char name[NAME_MAX + 1];
  pid_t owner;
  CUipcMemHandle handle;
  sem_t var_mutex;
} sh_dptr;

typedef struct sh_cu_event_st {
  char name[NAME_MAX + 1];
  pid_t owner;
  CUipcEventHandle handle;
  sem_t var_mutex;
} sh_event;

typedef struct sh_buffer_st {
  char name[NAME_MAX + 1];
  bool isBoard;
  size_t size;
  size_t data_size;
  uint32_t nb_proc;
  void *p_shm;
  sem_t mutex;
  sem_t wait_pub;
} sh_buffer;

typedef struct sh_barrier_st {
  char name[NAME_MAX + 1];
  bool valid;
  uint32_t nb_proc;
  uint32_t waiters_cnt;
  uint32_t val;
  sem_t b_sem;
  sem_t var_mutex;
} sh_barrier;

class CarmaIPCS {
 private:
  /*
     maps used to control the shms
   */
  std::map<uint32_t, sh_dptr *> dptrs;
  std::map<uint32_t, sh_event *> events;
  std::map<uint32_t, sh_barrier *> barriers;
  std::map<uint32_t, sh_buffer *> buffers;

  /*
     General purpose methods
   */
  void *create_shm(const char *name, size_t size);
  void *get_shm(const char *name);
  void free_shm(const char *name, void *p_shm, size_t size);
  void close_shm(void *p_shm, size_t size);
  void complete_clean();

  /*
    Transfer via CPU memory private methods
  */
  sh_buffer *get_elem_tshm(uint32_t id);
  int32_t write_gpu(void *dst, CUdeviceptr src, size_t bsize);
  int32_t read_gpu(CUdeviceptr dst, void *src, size_t bsize);

 public:
  /*
   general purpose methods
  */
  CarmaIPCS();
  ~CarmaIPCS();

  /*
    Cuda handles methods
  */
  // register a device pointer, id must be a non nul argument
  int32_t register_cudptr(uint32_t id, CUdeviceptr dptr);
  // register an event, id must be a non nul argument
  int32_t register_cuevent(uint32_t id, CUevent event);

  // get a memory handle
  int32_t get_memHandle(uint32_t id, CUipcMemHandle *phandle);
  // get a event handle
  int32_t get_eventHandle(uint32_t id, CUipcEventHandle *phandle);

  // free a memory handle shared mem space
  void free_memHandle(uint32_t id);
  // free a event handle shared event space
  void free_eventHandle(uint32_t id);

  /*
    Transfer via CPU memory methods
  */
  // allocation of the shm for memory tranfers
  int32_t alloc_transfer_shm(uint32_t id, size_t bsize, bool isBoard = false);
  // return size in bytes of the transfer shm in bsize
  int32_t get_size_transfer_shm(uint32_t id, size_t *bsize);
  // return actual data size used in bytes of the transfer shm in bsize
  int32_t get_datasize_transfer_shm(uint32_t id, size_t *bsize);
  // map the shm buffer to the cuda device
  int32_t map_transfer_shm(uint32_t id);
  // write to a transfer shm referenced by id
  int32_t write_transfer_shm(uint32_t id, const void *src, size_t bsize,
                         bool gpuBuffer = false);
  // reads from a transfer shm referenced by id
  int32_t read_transfer_shm(uint32_t id, void *dst, size_t bsize,
                        bool gpuBuffer = false);
  // map the shm buffer to the cuda device
  int32_t unmap_transfer_shm(uint32_t id);
  // free transfer shm ref by id
  void free_transfer_shm(uint32_t id);

  /*
    Barrier methods
  */
  // initialize a barrier with value the number of process who will call
  // wait_barrier
  int32_t init_barrier(uint32_t id, uint32_t value);
  // wait for the other process subscribed to the barrier
  int32_t wait_barrier(uint32_t id);
  // free the barrier structure, all blocked process will be unlocked<
  void free_barrier(uint32_t id);
};

#endif  // CARMA_IPCS_H
