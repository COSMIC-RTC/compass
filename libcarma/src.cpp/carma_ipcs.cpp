// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
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

//! \file      carma_ipc.cpp
//! \ingroup   libcarma
//! \class     carma_ipc
//! \brief     this class provides the ipc features to CarmaObj
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.3.0
//! \date      2022/01/24
//! \copyright GNU Lesser General Public License

#if 0
#include <carma_ipcs.h>

/*
 private methods
*/
void * CarmaIPCS::create_shm(const char *name, size_t size) {
  int tmp_fd;
  void * ret_ptr;

  if((tmp_fd = shm_open(name, O_CREAT | O_EXCL | O_RDWR, 0600)) == -1) //shm_open error, errno already set
    return NULL;
  if(ftruncate(tmp_fd, size) == -1) { //ftruncate error, same
    close(tmp_fd);
    return NULL;
  }
  ret_ptr = mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED, tmp_fd, 0);
  //if error, return value is consistent and errno is set
  close(tmp_fd);
  return ret_ptr;
}


void * CarmaIPCS::get_shm(const char *name) {
  int tmp_fd;
  void * ret_ptr;
  struct stat buf;

  if((tmp_fd = shm_open(name, O_RDWR, 0)) == -1) //shm_open error, errno already set
    return NULL;
  if(fstat(tmp_fd, &buf) == -1) {
    close(tmp_fd);
    return NULL;
  }
  ret_ptr = mmap(NULL, buf.st_size, PROT_READ | PROT_WRITE, MAP_SHARED, tmp_fd, 0);
  //if error, return value is consistent and errno is set
  close(tmp_fd);
  return ret_ptr;
}

void CarmaIPCS::free_shm(const char *name, void *p_shm, size_t size) {
  shm_unlink(name);
  munmap(p_shm, size);
}

void CarmaIPCS::close_shm(void *p_shm, size_t size) {
  munmap(p_shm, size);
}


void CarmaIPCS::complete_clean() {
  if(!dptrs.empty()) {
    std::map<unsigned int, sh_dptr *>::iterator it;
    for(it = dptrs.begin(); it!=dptrs.end(); ++it)
      free_memHandle(it->first);
    dptrs.clear();
  }
  if(!events.empty()) {
    std::map<unsigned int, sh_event *>::iterator it;
    for(it = events.begin(); it!=events.end(); ++it)
      free_eventHandle(it->first);
    events.clear();
  }
  if(!buffers.empty()) {
    std::map<unsigned int, sh_buffer *>::iterator it;
    for(it = buffers.begin(); it!=buffers.end(); ++it)
      free_transfer_shm(it->first);
    buffers.clear();
  }
  if(!barriers.empty()) {
    std::map<unsigned int, sh_barrier *>::iterator it;
    for(it = barriers.begin(); it!=barriers.end(); ++it)
      free_barrier(it->first);
    barriers.clear();
  }
}



/*
 general purpose methods
*/
CarmaIPCS::CarmaIPCS():
  dptrs(), events(), barriers(), buffers() {

}

CarmaIPCS::~CarmaIPCS() {
  complete_clean();
}



/*
  Cuda handles methods
*/

int CarmaIPCS::register_cudptr(unsigned int id,CUdeviceptr _dptr) {
  int res = EXIT_FAILURE;
  CUresult err;

  std::map<unsigned int, sh_dptr *>::iterator it;
  it = dptrs.find(id);
  if(it != dptrs.end()) { //id found in map
    errno = EEXIST;
    return res;
  }

  //create a new slot
  sh_dptr *dptr;
  char name[NAME_MAX+1];

  sprintf(name, "/cipcs_cudptr_%d", id);
  if((dptr = (sh_dptr *) create_shm(name, sizeof(sh_dptr))) == NULL)
    return res;
  //shm ok, structure initialization
  err = cuIpcGetMemHandle(&dptr->handle, _dptr);
  if(err != CUDA_SUCCESS) {
    free_shm(name, dptr, sizeof(sh_dptr));
    errno = err;
    res = -2; //to handle CUDA errors, see CIPCS_CHECK macro
    return res;
  }
  if(sem_init(&dptr->var_mutex, 1, 0) == -1) {
    free_shm(name, dptr, sizeof(sh_dptr));
    return res;
  }

  dptr->owner = getpid();
  strcpy(dptr->name, name);

  sem_post(&dptr->var_mutex); //mutex = 1;

  dptrs[id] = dptr;
  errno = 0;
  res = EXIT_SUCCESS;
  return res;
}


int CarmaIPCS::register_cuevent(unsigned int id, CUevent _event) {
  int res = EXIT_FAILURE;
  CUresult err;

  std::map<unsigned int, sh_event *>::iterator it;
  it = events.find(id);
  if(it != events.end()) { //id found in map
    errno = EEXIST;
    return res;
  }

  //create a new slot
  sh_event *event;
  char name[NAME_MAX+1];

  sprintf(name, "/cipcs_cuevent_%d", id);
  if((event = (sh_event *) create_shm(name, sizeof(sh_event))) == NULL)
    return res;
  //shm ok, structure initialization
  err = cuIpcGetEventHandle(&event->handle, _event);
  if(err != CUDA_SUCCESS) {
    free_shm(name, event, sizeof(sh_event));
    errno = err;
    res = -2; //to handle CUDA errors, see CIPCS_CHECK macro
    return res;
  }
  if(sem_init(&event->var_mutex, 1, 0) == -1) {
    free_shm(name, event, sizeof(sh_event));
    return res;
  }

  event->owner = getpid();
  strcpy(event->name, name);

  sem_post(&event->var_mutex); //mutex = 1;

  events[id] = event;

  errno = 0;
  res = EXIT_SUCCESS;
  return res;
}




int CarmaIPCS::get_memHandle(unsigned int id, CUipcMemHandle *phandle) {
  int res = EXIT_FAILURE;

  std::map<unsigned int, sh_dptr *>::iterator it;
  it = dptrs.find(id);
  if(it == dptrs.end()) { //not found in map
    sh_dptr *dptr;
    char name[NAME_MAX+1];
    sprintf(name, "/cipcs_cudptr_%d", id);
    if((dptr = (sh_dptr *)get_shm(name)) == NULL) { //not found in shm, id invalid
      errno = EINVAL;
      return res;
    } else { //found in shm
      if(sem_wait(&dptr->var_mutex))//probably destroyed while waiting on it, i.e. no more sharing
        return res;
      if(dptr->owner == getpid()) { //owner try to get handler, illegal operation
        sem_post(&dptr->var_mutex);
        errno = EPERM;
        return res;
      }
      //else, add it to the local map
      sem_post(&dptr->var_mutex);
      dptrs[id] = dptr;
    }
  }
  //everything's OK
  *phandle = dptrs[id]->handle;

  errno = 0;
  res = EXIT_SUCCESS;
  return res;
}



int CarmaIPCS::get_eventHandle(unsigned int id, CUipcEventHandle *phandle) {
  int res = EXIT_FAILURE;

  std::map<unsigned int, sh_event *>::iterator it;
  it = events.find(id);
  if(it == events.end()) { //not found in map
    sh_event *event;
    char name[NAME_MAX+1];
    sprintf(name, "/cipcs_cuevent_%d", id);
    if((event = (sh_event *)get_shm(name)) == NULL) { //not found in shm
      errno = EINVAL;
      return res;
    } else { //found in shm
      if(sem_wait(&event->var_mutex))//probably destroyed while waiting on it, i.e. no more sharing
        return res;
      if(event->owner == getpid()) { //owner try to get handler, illegal operation
        sem_post(&event->var_mutex);
        errno = EPERM;
        return res;
      }
      //else, add it to the local map
      sem_post(&event->var_mutex);
      events[id] = event;
    }
  }
  //else, everything's OK
  *phandle = events[id]->handle;

  errno = 0;
  res = EXIT_SUCCESS;
  return res;
}


void CarmaIPCS::free_memHandle(unsigned int id) {
  std::map<unsigned int, sh_dptr *>::iterator it;
  it = dptrs.find(id);
  if(it == dptrs.end()) { //not found in map
    return;
  }
  sh_dptr *dptr = dptrs[id];
  dptrs.erase(it);
  sem_wait(&dptr->var_mutex);
  if(dptr->owner == getpid()) {
    sem_destroy(&dptr->var_mutex);
    free_shm(dptr->name, dptr, sizeof(sh_dptr));
  } else {
    sem_post(&dptr->var_mutex);
    close_shm(dptr, sizeof(sh_dptr));
  }
}


void CarmaIPCS::free_eventHandle(unsigned int id) {
  std::map<unsigned int, sh_event *>::iterator it;
  it = events.find(id);
  if(it == events.end()) { //not found in map
    return;
  }
  sh_event *event = events[id];
  events.erase(it);
  sem_wait(&event->var_mutex);
  if(event->owner == getpid()) {
    sem_destroy(&event->var_mutex);
    free_shm(event->name, event, sizeof(sh_event));
  } else {
    sem_post(&event->var_mutex);
    close_shm(event, sizeof(sh_event));
  }
}




/*
  Transfer via CPU memory methods
*/
sh_buffer *CarmaIPCS::get_elem_tshm(unsigned int id) {
  std::map<unsigned int, sh_buffer *>::iterator it;

  it = buffers.find(id);
  if(it == buffers.end()) { //not found in map
    //STAMP("tranfer not found in map\n");
    char name[NAME_MAX+1];
    sh_buffer *buffer;
    sprintf(name, "/cipcs_cubuffer_%d", id);
    if((buffer = (sh_buffer *)get_shm(name)) == NULL) { //not found in shm
      errno = EINVAL;
      return NULL;
    } else { //found in shm
      //get included shm
      //STAMP("tranfer found in shm\n");
      char in_name[NAME_MAX+1];
      sprintf(in_name, "/cipcs_cubuffer_%d_in", id);
      if((buffer->p_shm = get_shm(in_name)) == NULL) {
        close_shm(buffer, sizeof(sh_buffer));
        return NULL; //errno set by get_shm
      }

      if(sem_wait(&(buffer->mutex)) == -1)//probably destroyed while waiting on it, i.e. no more sharing
        return NULL;
      //else, add it to the local map
      buffer->nb_proc++;
      sem_post(&(buffer->mutex));
      buffers[id] = buffer;
    }
  }

  return buffers[id];
}



int CarmaIPCS::alloc_transfer_shm(unsigned int id, size_t bsize, bool isBoard) {
  int res = EXIT_FAILURE;
  std::map<unsigned int, sh_buffer *>::iterator it;
  it = buffers.find(id);
  if(it != buffers.end()) { //id found in map
    errno = EEXIST;
    return res;
  }

  //create a new slot
  sh_buffer *buffer=NULL;
  char name[NAME_MAX+1];
  char in_name[NAME_MAX+1];

  sprintf(name, "/cipcs_cubuffer_%d", id);
  sprintf(in_name, "/cipcs_cubuffer_%d_in", id);
  if((buffer = (sh_buffer *) create_shm(name, sizeof(sh_buffer))) == NULL)
    return res;
  //shm ok, struct initialization
  if((buffer->p_shm = create_shm(in_name, bsize)) == NULL) {
    free_shm(name, buffer, sizeof(sh_buffer));
    return res;
  }
  if(sem_init(&(buffer->mutex), 1, 0) == -1) {
    free_shm(in_name, buffer->p_shm, bsize);
    free_shm(name, buffer, sizeof(sh_buffer));
    return res;
  }
  if(isBoard)
    if(sem_init(&(buffer->wait_pub), 1, 0) == -1) {
      sem_destroy(&(buffer->mutex));
      free_shm(in_name, buffer->p_shm, bsize);
      free_shm(name, buffer, sizeof(sh_buffer));
      return res;
    }

  buffer->size = bsize;
  buffer->data_size = 0;
  buffer->nb_proc = 1;
  strcpy(buffer->name, name);
  buffer->isBoard = isBoard;

  sem_post(&(buffer->mutex)); //mutex = 1;

  buffers[id] = buffer;

  errno = 0;
  res = EXIT_SUCCESS;
  return res;
}


int CarmaIPCS::get_size_transfer_shm(unsigned int id, size_t *bsize) {
  int res = EXIT_FAILURE;
  sh_buffer *buffer=NULL;

  if( (buffer = get_elem_tshm(id)) == NULL)
    return res; //errno previously set

  *bsize = buffer->size;

  errno = 0;
  res = EXIT_SUCCESS;
  return res;
}


int CarmaIPCS::get_datasize_transfer_shm(unsigned int id, size_t *bsize) {
  int res = EXIT_FAILURE;
  sh_buffer *buffer=NULL;
  //STAMP("in get_datasize\n");
  if( (buffer = get_elem_tshm(id)) == NULL)
    return res; //errno previously set

  //STAMP("\tsem_wait, mutex\n");
  if(sem_wait(&(buffer->mutex)) == -1)
    return res;
  *bsize = buffer->data_size;

  //STAMP("\t sem_post, mutex\n");
  sem_post(&(buffer->mutex));

  errno = 0;
  res = EXIT_SUCCESS;
  //STAMP("out get_datasize\n");
  return res;
}



int CarmaIPCS::map_transfer_shm(unsigned int id) {
  int res = EXIT_FAILURE;
  sh_buffer *buffer=NULL;

  if( (buffer = get_elem_tshm(id)) == NULL)
    return res; //errno previously set

  cuMemHostRegister(buffer->p_shm, buffer->data_size, CU_MEMHOSTREGISTER_DEVICEMAP);

  errno = 0;
  res = EXIT_SUCCESS;
  return res;
}



int CarmaIPCS::write_gpu(void *dst, CUdeviceptr src, size_t bsize) {
  int res = EXIT_FAILURE;

  if((errno = cuMemcpyDtoH(dst, src, bsize)) != CUDA_SUCCESS)
    res = -2;
  else {
    errno = 0;
    res = EXIT_SUCCESS;
  }
  return res;
}


int CarmaIPCS::read_gpu(CUdeviceptr dst, void *src, size_t bsize) {
  int res = EXIT_FAILURE;

  if((errno = cuMemcpyHtoD(dst, src, bsize)) != CUDA_SUCCESS)
    res = -2;
  else {
    errno = 0;
    res = EXIT_SUCCESS;
  }
  return res;
}


int CarmaIPCS::write_transfer_shm(unsigned int id, const void *src, size_t bsize, bool isGpuBuffer) {
  int res = EXIT_FAILURE;
  sh_buffer *buffer=NULL;

  //STAMP("in write\n");
  if( (buffer = get_elem_tshm(id)) == NULL)
    return res; //errno previously set

  //STAMP("\tsem_wait, mutex\n");
  if(sem_wait(&(buffer->mutex)) == -1)
    return res;

  //todo : error code, maybe?
  if(bsize > buffer->size)
    bsize = buffer->size;

  if(isGpuBuffer) {
    //STAMP("\t\tdevice buffer\n");
    if((res = write_gpu(buffer->p_shm, (CUdeviceptr) src, bsize)) < 0) {
      sem_post(&(buffer->mutex));
      return res;
    }
    //STAMP("\t\tCopy ok\n");
  } else {
    //STAMP("\t\thost buffer\n");
    memcpy(buffer->p_shm, src, bsize);
    //STAMP("\t\tCopy ok\n");
  }

  buffer->data_size = bsize;

  if(buffer->isBoard)
    sem_post(&(buffer->wait_pub));

  //STAMP("\tsem_post mutex\n");
  sem_post(&(buffer->mutex));
  errno = 0;
  res = EXIT_SUCCESS;
  //STAMP("out write\n");
  return res;
}


int CarmaIPCS::read_transfer_shm(unsigned int id, void *dst, size_t bsize, bool isGpuBuffer) {
  int res = EXIT_FAILURE;
  sh_buffer *buffer=NULL;

  //STAMP("in read\n");

  if( (buffer = get_elem_tshm(id)) == NULL)
    return res; //errno previously set

  //STAMP("\tsem_wait, mutex\n");
  if(sem_wait(&(buffer->mutex)) == -1)
    return res;

  //todo : error code, maybe?
  if(bsize > buffer->size)
    bsize = buffer->size;

  if(buffer->isBoard) {
    sem_post(&(buffer->mutex));
    sem_wait(&(buffer->wait_pub));
    //todo: to improve, this could be buggy.
    sem_wait(&(buffer->mutex));
  }
  if(isGpuBuffer) {

    //STAMP("\t\t device\n");
    if((res = read_gpu((CUdeviceptr) dst, buffer->p_shm, bsize)) < 0) {
      sem_post(&(buffer->mutex));
      return res;
    }
  } else {
    //STAMP("\t\t host\n");
    memcpy(dst, buffer->p_shm, bsize);
  }


  //STAMP("\tsem_post, mutex\n");
  sem_post(&(buffer->mutex));
  errno = 0;
  res = EXIT_SUCCESS;

  //STAMP("out read\n");
  return res;
}

int CarmaIPCS::unmap_transfer_shm(unsigned int id) {
  int res = EXIT_FAILURE;
  sh_buffer *buffer=NULL;

  if( (buffer = get_elem_tshm(id)) == NULL)
    return res; //errno previously set

  cuMemHostUnregister(buffer->p_shm);

  errno = 0;
  res = EXIT_SUCCESS;
  return res;
}


void CarmaIPCS::free_transfer_shm(unsigned int id) {
  std::map<unsigned int, sh_buffer *>::iterator it;

  it = buffers.find(id);
  if(it == buffers.end()) { //not found in map
    return;
  }
  sh_buffer *buffer = buffers[id];
  buffers.erase(it);
  sem_wait(&(buffer->mutex));
  buffer->nb_proc--;
  if(buffer->nb_proc == 0) {
    if(buffer->isBoard)
      sem_destroy(&(buffer->wait_pub));
    sem_destroy(&(buffer->mutex));
    char in_name[NAME_MAX+4];
    sprintf(in_name, "%s_in", buffer->name);

    free_shm(in_name, buffer->p_shm, buffer->size);
    free_shm(buffer->name, buffer, sizeof(sh_buffer));
  } else {
    sem_post(&(buffer->mutex));
    close_shm(buffer->p_shm, buffer->size);
    close_shm(buffer, sizeof(sh_buffer));
  }
}




/*
  Barrier methods
*/
int CarmaIPCS::init_barrier(unsigned int id, unsigned int value) {
  int res = EXIT_FAILURE;
  std::map<unsigned int, sh_barrier *>::iterator it;
  it = barriers.find(id);
  if(it != barriers.end()) { //id found in map
    errno = EEXIST;
    return res;
  }

  //create a new slot
  sh_barrier *barrier;
  char name[NAME_MAX+1];

  sprintf(name, "/cipcs_cubarrier_%d", id);
  if((barrier = (sh_barrier *) create_shm(name, sizeof(sh_barrier))) == NULL)
    return res;
  //shm ok, struct initialization
  if(sem_init(&barrier->var_mutex, 1, 0) == -1) {
    free_shm(name, barrier, sizeof(sh_barrier));
    return res;
  }
  if(sem_init(&barrier->b_sem, 1, 0) == -1) {
    sem_destroy(&barrier->var_mutex);
    free_shm(name, barrier, sizeof(sh_barrier));
    return res;
  }

  barrier->valid = true;
  barrier->waiters_cnt = 0;
  barrier->nb_proc = 1;
  barrier->val = value;
  strcpy(barrier->name, name);

  sem_post(&barrier->var_mutex); //mutex = 1;

  barriers[id] = barrier;

  errno = 0;
  res = EXIT_SUCCESS;
  return res;
}



int CarmaIPCS::wait_barrier(unsigned int id) {
  int res = EXIT_FAILURE;
  std::map<unsigned int, sh_barrier *>::iterator it;
  sh_barrier *barrier;

  it = barriers.find(id);
  if(it == barriers.end()) { //not found in map
    //STAMP("%d, barrier not found in map\n", getpid());
    char name[NAME_MAX+1];

    sprintf(name, "/cipcs_cubarrier_%d", id);
    if((barrier = (sh_barrier *)get_shm(name)) == NULL) { //not found in shm
      //STAMP("%d, barrier not found in shm\n", getpid());
      errno = EINVAL;
      return res;
    }
    //else, found in shm, add to map
    //STAMP("%d, wait for var_mutex 1 (adding)\n", getpid());
    sem_wait(&barrier->var_mutex);
    barrier->nb_proc++;
    barriers[id] = barrier;
    sem_post(&barrier->var_mutex);
    //STAMP("%d, barrier added to map\n", getpid());
  } else
    barrier = barriers[id];

  //access barrier
  if(sem_wait(&barrier->var_mutex) == -1)
    return res;
  if(!barrier->valid) {
    //STAMP("%d, post var_mutex invalid\n", getpid());
    sem_post(&barrier->var_mutex);
    errno = EBADF;
    return res;
  }
  barrier->waiters_cnt++;
//STAMP("waiters_cnt %d, val %d\n", barrier->waiters_cnt, barrier->val);
  if(barrier->waiters_cnt == barrier->val) {
    //STAMP("%d, post var_mutex, unlock\n", getpid());
    sem_post(&barrier->var_mutex);
    for(unsigned int i = 1; i < barrier->val; ++i) {
      sem_post(&barrier->b_sem); //(val - 1) post
    }
  } else { //release var mutex and wait on barrier semaphore
    //STAMP("%d, post var_mutex, wait b_sem\n", getpid());
    sem_post(&barrier->var_mutex);
    sem_wait(&barrier->b_sem);
    //STAMP("%d, b_sem unlocked\n", getpid());
    if(!barrier->valid) {
      //STAMP("%d, Barrier invalidated\n", getpid());
      errno = ECANCELED;
      return res;
    }
  }
  sem_wait(&barrier->var_mutex);
  barrier->waiters_cnt--;
  sem_post(&barrier->var_mutex);

  return res;
}


void CarmaIPCS::free_barrier(unsigned int id) {
  std::map<unsigned int, sh_barrier *>::iterator it;
  it = barriers.find(id);
  if(it == barriers.end()) { //not found in map
    return;
  }
  sh_barrier *barrier = barriers[id];
  barriers.erase(it);
  //STAMP("%d, wait for var_mutex\n", getpid());
  sem_wait(&barrier->var_mutex);
  barrier->nb_proc--;
  barrier->valid = false;
  if(barrier->waiters_cnt) {
    for(unsigned int i = 0; i < barrier->waiters_cnt; ++i)
      sem_post(&barrier->b_sem);
  }
  if(barrier->nb_proc == 0) {
    sem_destroy(&barrier->var_mutex);
    sem_destroy(&barrier->b_sem);
    free_shm(barrier->name, barrier, sizeof(sh_barrier));
  } else {
    sem_post(&barrier->var_mutex);
    close_shm(barrier, sizeof(sh_barrier));
  }
}
#endif
