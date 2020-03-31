// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2019 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.
//  Distributed under GNU - LGPL
//
//  COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
//  General Public License as published by the Free Software Foundation, either version 3 of the License,
//  or any later version.
//
//  COMPASS: End-to-end AO simulation tool using GPU acceleration
//  The COMPASS platform was designed to meet the need of high-performance for the simulation of AO systems.
//
//  The final product includes a software package for simulating all the critical subcomponents of AO,
//  particularly in the context of the ELT and a real-time core based on several control approaches,
//  with performances consistent with its integration into an instrument. Taking advantage of the specific
//  hardware architecture of the GPU, the COMPASS tool allows to achieve adequate execution speeds to
//  conduct large simulation campaigns called to the ELT.
//
//  The COMPASS platform can be used to carry a wide variety of simulations to both testspecific components
//  of AO of the E-ELT (such as wavefront analysis device with a pyramid or elongated Laser star), and
//  various systems configurations such as multi-conjugate AO.
//
//  COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
//  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//  See the GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License along with COMPASS.
//  If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>.
// -----------------------------------------------------------------------------

//! \file      carma_ipc.h
//! \ingroup   libcarma
//! \class     carma_ipc
//! \brief     this class provides the ipc features to carma_obj
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.4.1
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

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
#include <sys/types.h>
#include <unistd.h>

#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>

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
  unsigned int nb_proc;
  void *p_shm;
  sem_t mutex;
  sem_t wait_pub;
} sh_buffer;

typedef struct sh_barrier_st {
  char name[NAME_MAX + 1];
  bool valid;
  unsigned int nb_proc;
  unsigned int waiters_cnt;
  unsigned int val;
  sem_t b_sem;
  sem_t var_mutex;
} sh_barrier;

class carma_ipcs {
 private:
  /*
     maps used to control the shms
   */
  std::map<unsigned int, sh_dptr *> dptrs;
  std::map<unsigned int, sh_event *> events;
  std::map<unsigned int, sh_barrier *> barriers;
  std::map<unsigned int, sh_buffer *> buffers;

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
  sh_buffer *get_elem_tshm(unsigned int id);
  int write_gpu(void *dst, CUdeviceptr src, size_t bsize);
  int read_gpu(CUdeviceptr dst, void *src, size_t bsize);

 public:
  /*
   general purpose methods
  */
  carma_ipcs();
  ~carma_ipcs();

  /*
    Cuda handles methods
  */
  // register a device pointer, id must be a non nul argument
  int register_cudptr(unsigned int id, CUdeviceptr dptr);
  // register an event, id must be a non nul argument
  int register_cuevent(unsigned int id, CUevent event);

  // get a memory handle
  int get_memHandle(unsigned int id, CUipcMemHandle *phandle);
  // get a event handle
  int get_eventHandle(unsigned int id, CUipcEventHandle *phandle);

  // free a memory handle shared mem space
  void free_memHandle(unsigned int id);
  // free a event handle shared event space
  void free_eventHandle(unsigned int id);

  /*
    Transfer via CPU memory methods
  */
  // allocation of the shm for memory tranfers
  int alloc_transfer_shm(unsigned int id, size_t bsize, bool isBoard = false);
  // return size in bytes of the transfer shm in bsize
  int get_size_transfer_shm(unsigned int id, size_t *bsize);
  // return actual data size used in bytes of the transfer shm in bsize
  int get_datasize_transfer_shm(unsigned int id, size_t *bsize);
  // map the shm buffer to the cuda device
  int map_transfer_shm(unsigned int id);
  // write to a transfer shm referenced by id
  int write_transfer_shm(unsigned int id, const void *src, size_t bsize,
                         bool gpuBuffer = false);
  // reads from a transfer shm referenced by id
  int read_transfer_shm(unsigned int id, void *dst, size_t bsize,
                        bool gpuBuffer = false);
  // map the shm buffer to the cuda device
  int unmap_transfer_shm(unsigned int id);
  // free transfer shm ref by id
  void free_transfer_shm(unsigned int id);

  /*
    Barrier methods
  */
  // initialize a barrier with value the number of process who will call
  // wait_barrier
  int init_barrier(unsigned int id, unsigned int value);
  // wait for the other process subscribed to the barrier
  int wait_barrier(unsigned int id);
  // free the barrier structure, all blocked process will be unlocked<
  void free_barrier(unsigned int id);
};

#endif  // CARMA_IPCS_H
