#include <carma_ipcs.h>

/*
 private methods
*/
void * carma_ipcs::create_shm(const char *name, size_t size){
  int * res = (int *)EXIT_FAILURE;
  int tmp_fd;
  void * ret_ptr;

  if((tmp_fd = shm_open(name, O_CREAT | O_EXCL | O_RDWR, 0600)) == -1) //shm_open error, errno already set
    return res;
  if(ftruncate(tmp_fd, size) == -1){ //ftruncate error, same
    close(tmp_fd);
    return res;
  }
  ret_ptr = mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED, tmp_fd, 0);
  //if error, return value is consistent and errno is set
  close(tmp_fd);
  return ret_ptr;
}


void * carma_ipcs::get_shm(const char *name){
  int * res = (int *) EXIT_FAILURE;
  int tmp_fd;
  void * ret_ptr;
  struct stat buf;

  if((tmp_fd = shm_open(name, O_RDWR, 0)) == -1) //shm_open error, errno already set
    return res;
  if(fstat(tmp_fd, &buf) == -1){
    close(tmp_fd);
    return res;
  }
  ret_ptr = mmap(NULL, buf.st_size, PROT_READ | PROT_WRITE, MAP_SHARED, tmp_fd, 0);
  //if error, return value is consistent and errno is set
  close(tmp_fd);
  return ret_ptr;
}

void carma_ipcs::free_shm(const char *name, void *p_shm, size_t size){
  shm_unlink(name);
  munmap(p_shm, size);
}

void carma_ipcs::close_shm(void *p_shm, size_t size){
  munmap(p_shm, size);
}


void carma_ipcs::complete_clean(){
  if(!dptrs.empty()){
    std::map<unsigned int, sh_dptr *>::iterator it;
    for(it = dptrs.begin(); it!=dptrs.end(); ++it)
      free_memHandle(it->first);
    dptrs.clear();
  }
  if(!events.empty()){
    std::map<unsigned int, sh_event *>::iterator it;
    for(it = events.begin(); it!=events.end(); ++it)
      free_eventHandle(it->first);
    events.clear();
  }
  if(!buffers.empty()){
    std::map<unsigned int, sh_buffer *>::iterator it;
    for(it = buffers.begin(); it!=buffers.end(); ++it)
      free_transfer_shm(it->first);
    buffers.clear();
  }
  if(!barriers.empty()){
    std::map<unsigned int, sh_barrier *>::iterator it;
    for(it = barriers.begin(); it!=barriers.end(); ++it)
      free_barrier(it->first);
    barriers.clear();
  }
}



/*
 general purpose methods
*/
carma_ipcs::carma_ipcs():
dptrs(), events(), barriers(), buffers()
{

}

carma_ipcs::~carma_ipcs(){
  complete_clean();
}



/*
  Cuda handles methods
*/

int carma_ipcs::register_cudptr(unsigned int id ,CUdeviceptr _dptr){
  int res = EXIT_FAILURE;
  CUresult err;

  std::map<unsigned int, sh_dptr *>::iterator it;
  it = dptrs.find(id);
  if(it != dptrs.end()){ //id found in map
    errno = EEXIST;
    return res;
  }

  //create a new slot
  sh_dptr *dptr;
  char name[NAME_MAX+1];

  sprintf(name, "/cipcs_cudptr_%d", id);
  if((dptr = (sh_dptr *) create_shm(name, sizeof(sh_dptr))) == (void *) -1)
    return res;
  //shm ok, structure initialization
  err = cuIpcGetMemHandle(&dptr->handle, _dptr);
  if(err != CUDA_SUCCESS){
    free_shm(name, dptr, sizeof(sh_dptr));
    errno = err;
    res = -2; //to handle CUDA errors, see CIPCS_CHECK macro
    return res;
  }
  if(sem_init(&dptr->var_mutex, 1, 0) == -1){
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


int carma_ipcs::register_cuevent(unsigned int id, CUevent _event){
  int res = EXIT_FAILURE;
  CUresult err;

  std::map<unsigned int, sh_event *>::iterator it;
  it = events.find(id);
  if(it != events.end()){ //id found in map
    errno = EEXIST;
    return res;
  }

  //create a new slot
  sh_event *event;
  char name[NAME_MAX+1];

  sprintf(name, "/cipcs_cuevent_%d", id);
  if((event = (sh_event *) create_shm(name, sizeof(sh_event))) == (void *) -1)
    return res;
  //shm ok, structure initialization
  err = cuIpcGetEventHandle(&event->handle, _event);
  if(err != CUDA_SUCCESS){
    free_shm(name, event, sizeof(sh_event));
    errno = err;
    res = -2; //to handle CUDA errors, see CIPCS_CHECK macro
    return res;
  }
  if(sem_init(&event->var_mutex, 1, 0) == -1){
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




int carma_ipcs::get_memHandle(unsigned int id, CUipcMemHandle *phandle){
  int res = EXIT_FAILURE;

  std::map<unsigned int, sh_dptr *>::iterator it;
  it = dptrs.find(id);
  if(it == dptrs.end()){ //not found in map
    sh_dptr *dptr;
    char name[NAME_MAX+1];
    sprintf(name, "/cipcs_cudptr_%d", id);
    if((dptr = (sh_dptr *)get_shm(name)) == (void *) -1){ //not found in shm, id invalid
      errno = EINVAL;
      return res;
    }
    else{ //found in shm
      if(sem_wait(&dptr->var_mutex))//probably destroyed while waiting on it, i.e. no more sharing
        return res;
      if(dptr->owner == getpid()){ //owner try to get handler, illegal operation
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



int carma_ipcs::get_eventHandle(unsigned int id, CUipcEventHandle *phandle){
  int res = EXIT_FAILURE;

  std::map<unsigned int, sh_event *>::iterator it;
  it = events.find(id);
  if(it == events.end()){ //not found in map
    sh_event *event;
    char name[NAME_MAX+1];
    sprintf(name, "/cipcs_cuevent_%d", id);
    if((event = (sh_event *)get_shm(name)) == (void *) -1){ //not found in shm
      errno = EINVAL;
      return res;
    }
    else{ //found in shm
      if(sem_wait(&event->var_mutex))//probably destroyed while waiting on it, i.e. no more sharing
        return res;
      if(event->owner == getpid()){ //owner try to get handler, illegal operation
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


void carma_ipcs::free_memHandle(unsigned int id){
  std::map<unsigned int, sh_dptr *>::iterator it;
  it = dptrs.find(id);
  if(it == dptrs.end()){ //not found in map
    return;
  }
  sh_dptr *dptr = dptrs[id];
  dptrs.erase(it);
  sem_wait(&dptr->var_mutex);
  if(dptr->owner == getpid()){
    sem_destroy(&dptr->var_mutex);
    free_shm(dptr->name, dptr, sizeof(sh_dptr));
  }
  else{
    sem_post(&dptr->var_mutex);
    close_shm(dptr, sizeof(sh_dptr));
  }
}


void carma_ipcs::free_eventHandle(unsigned int id){
  std::map<unsigned int, sh_event *>::iterator it;
  it = events.find(id);
  if(it == events.end()){ //not found in map
    return;
  }
  sh_event *event = events[id];
  events.erase(it);
  sem_wait(&event->var_mutex);
  if(event->owner == getpid()){
    sem_destroy(&event->var_mutex);
    free_shm(event->name, event, sizeof(sh_event));
  }
  else{
    sem_post(&event->var_mutex);
    close_shm(event, sizeof(sh_event));
  }
}




/*
  Transfer via CPU memory methods
*/
int carma_ipcs::get_elem_tshm(unsigned int id, sh_buffer *buffer){
  int res = EXIT_FAILURE;
    std::map<unsigned int, sh_buffer *>::iterator it;

    it = buffers.find(id);
    if(it == buffers.end()){ //not found in map
      char name[NAME_MAX+1];
      sprintf(name, "/cipcs_cubuffer_%d", id);
      if((buffer = (sh_buffer *)get_shm(name)) == (void *) -1){ //not found in shm
        errno = EINVAL;
        return res;
      }
      else{ //found in shm
        //get included shm
        char name[NAME_MAX+1];
        sprintf(name, "/cipcs_cubuffer_%d_in", id);
        if((buffer->p_shm = get_shm(name)) == (void *) -1){
          close_shm(buffer, sizeof(sh_buffer));
          return res; //errno set by get_shm
        }

        if(sem_wait(&buffer->mutex) == -1)//probably destroyed while waiting on it, i.e. no more sharing
          return res;
        //else, add it to the local map
        buffer->nb_proc++;
        sem_post(&buffer->mutex);
        buffers[id] = buffer;
      }
    }
    else
      buffer = buffers[id];

    errno = 0;
    res = EXIT_SUCCESS;
    return res;
}



int carma_ipcs::alloc_transfer_shm(unsigned int id, size_t bsize, bool isBoard){
  int res = EXIT_FAILURE;
  std::map<unsigned int, sh_buffer *>::iterator it;
  it = buffers.find(id);
  if(it != buffers.end()){ //id found in map
    errno = EEXIST;
    return res;
  }

  //create a new slot
  sh_buffer *buffer;
  char name[NAME_MAX+1];
  char in_name[NAME_MAX+1];

  sprintf(name, "/cipcs_cubuffer_%d", id);
  sprintf(in_name, "/cipcs_cubuffer_%d_in", id);
  if((buffer = (sh_buffer *) create_shm(name, sizeof(sh_buffer))) == (void *) -1)
    return res;
  //shm ok, struct initialization
  if((buffer->p_shm = create_shm(in_name, bsize)) == (void *) -1){
    free_shm(name, buffer, sizeof(sh_buffer));
    return res;
  }
  if(sem_init(&buffer->mutex, 1, 0) == -1){
    free_shm(in_name, buffer->p_shm, bsize);
    free_shm(name, buffer, sizeof(sh_buffer));
    return res;
  }
  if(isBoard)
    if(sem_init(&buffer->wait_pub, 1, 0) == -1){
      sem_destroy(&buffer->mutex);
      free_shm(in_name, buffer->p_shm, bsize);
      free_shm(name, buffer, sizeof(sh_buffer));
      return res;
    }

  buffer->size = bsize;
  buffer->data_size = 0;
  buffer->nb_proc = 1;
  strcpy(buffer->name, name);
  buffer->isBoard = isBoard;

  sem_post(&buffer->mutex); //mutex = 1;

  buffers[id] = buffer;

  errno = 0;
  res = EXIT_SUCCESS;
  return res;
}


int get_size_transfer_shm(unsigned int id, size_t *bsize){
  int res = EXIT_FAILURE;
  sh_buffer *buffer = NULL;

  if(get_elem_tshm(id, buffer) == EXIT_FAILURE)
    return res; //errno previously set

  *bsize = buffer->size;

  errno = 0;
  res = EXIT_SUCCESS;
  return res;
}


int get_datasize_transfer_shm(unsigned int id, size_t *bsize){
  int res = EXIT_FAILURE;
  sh_buffer *buffer = NULL;

  if(get_elem_tshm(id, buffer) == EXIT_FAILURE)
    return res; //errno previously set

  if(sem_wait(&buffer->mutex) == -1)
    return res;
  *bsize = buffer->data_size;
  sem_post(&buffer->mutex);

  errno = 0;
  res = EXIT_SUCCESS;
  return res;
}



int carma_ipcs::map_transfer_shm(unsigned int id){
  int res = EXIT_FAILURE;
  sh_buffer *buffer = NULL;

  if(get_elem_tshm(id, buffer) == EXIT_FAILURE)
    return res; //errno previously set

  //TODO MAPPING CUDA

  errno = 0;
  res = EXIT_SUCCESS;
  return res;
}



int carma_ipcs::write_gpu(void *dst, CUdeviceptr src, size_t bsize){
  int res = EXIT_FAILURE;
  //todo
  errno = 0;
  res = EXIT_SUCCESS;
  return res;
}


int carma_ipcs::read_gpu(CUdeviceptr dst, void *src, size_t bsize){
  int res = EXIT_FAILURE;
  //todo
  errno = 0;
  res = EXIT_SUCCESS;
  return res;
}


int carma_ipcs::write_transfer_shm(unsigned int id, const void *src, size_t bsize, bool isGpuBuffer){
  int res = EXIT_FAILURE;
  sh_buffer *buffer = NULL;

  if(get_elem_tshm(id, buffer) == EXIT_FAILURE)
    return res; //errno previously set

  if(sem_wait(&buffer->mutex) == -1)
    return res;

  if(gpuBuffer){
    if((res = write_gpu(buffer->p_shm, (CUdeviceptr) src, bsize)) == -1){
      sem_post(&buffer->mutex);
      return res;
    }
  }
  else{
    memcpy(buffer->p_shm, src, bsize);
  }

  buffer->data_size = bsize;

  if(buffer->isBoard)
    sem_post(&buffer->wait_pub);

  sem_post(&buffer->mutex);
  errno = 0;
  res = EXIT_SUCCESS;
  return res;
}


int carma_ipcs::read_transfer_shm(unsigned int id, void *dst, size_t bsize, bool isGpuBuffer){
  int res = EXIT_FAILURE;
  sh_buffer *buffer = NULL;


  if(get_elem_tshm(id, buffer) == EXIT_FAILURE)
    return res; //errno previously set

  if(sem_wait(&buffer->mutex) == -1)
    return res;

  if(buffer->isBoard){
    sem_post(&buffer->mutex);
    sem_wait(&buffer->wait_pub);
  }

  if(gpuBuffer){
    if((res = read_gpu((CUdeviceptr) dst, buffer->p_shm, bsize)) == -1){
      sem_post(&buffer->mutex);
      return res;
    }
  }
  else{
    memcpy(src, buffer->p_shm, bsize);
  }

  sem_post(&buffer->mutex);
  errno = 0;
  res = EXIT_SUCCESS;
  return res;
}

int unmap_transfer_shm(unsigned int id){
  int res = EXIT_FAILURE;
  sh_buffer *buffer = NULL;

  if(get_elem_tshm(id, buffer) == EXIT_FAILURE)
    return res; //errno previously set

  //TODO MAPPING CUDA

  errno = 0;
  res = EXIT_SUCCESS;
  return res;
}


void carma_ipcs::free_transfer_shm(unsigned int id){
  std::map<unsigned int, sh_buffer *>::iterator it;

  it = buffers.find(id);
  if(it == buffers.end()){ //not found in map
    return;
  }
  sh_buffer *buffer = buffers[id];
  buffers.erase(it);
  sem_wait(&buffer->mutex);
  buffer->nb_proc--;
  if(buffer->nb_proc == 0){
    if(buffer->isBoard)
      sem_destroy(&buffer->wait_pub);
    sem_destroy(&buffer->mutex);
    char in_name[NAME_MAX+1];
    sprintf(in_name, "%s_in", buffer->name);

    free_shm(in_name, buffer->p_shm, buffer->size);
    free_shm(buffer->name, buffer, sizeof(sh_buffer));
  }
  else{
    sem_post(&buffer->mutex);
    close_shm(buffer->p_shm, buffer->size);
    close_shm(buffer, sizeof(sh_buffer));
  }
}




/*
  Barrier methods
*/
int carma_ipcs::init_barrier(unsigned int id, unsigned int value){
  int res = EXIT_FAILURE;

  std::map<unsigned int, sh_barrier *>::iterator it;
  it = barriers.find(id);
  if(it != barriers.end()){ //id found in map
    errno = EEXIST;
    return res;
  }

  //create a new slot
  sh_barrier *barrier;
  char name[NAME_MAX+1];

  sprintf(name, "/cipcs_cubarrier_%d", id);
  if((barrier = (sh_barrier *) create_shm(name, sizeof(sh_barrier))) == (void *) -1)
    return res;
  //shm ok, struct initialization
  if(sem_init(&barrier->var_mutex, 1, 0) == -1){
    free_shm(name, barrier, sizeof(sh_barrier));
    return res;
  }
  if(sem_init(&barrier->b_sem, 1, 0) == -1){
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



int carma_ipcs::wait_barrier(unsigned int id){
  int res = EXIT_FAILURE;
  std::map<unsigned int, sh_barrier *>::iterator it;
  sh_barrier *barrier;


  it = barriers.find(id);
  if(it == barriers.end()){ //not found in map
    //STAMP("%d, barrier not found in map\n", getpid());
    sh_barrier * barrier;
    char name[NAME_MAX+1];

    sprintf(name, "/cipcs_cubarrier_%d", id);
    if((barrier = (sh_barrier *)get_shm(name)) == (void *) -1){ //not found in shm
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
  }
  else
    barrier = barriers[id];

  //access barrier
  //STAMP("%d, wait for var_mutex 2 (access)\n", getpid());
  if(sem_wait(&barrier->var_mutex) == -1)
    return res;
  if(!barrier->valid){
    //STAMP("%d, post var_mutex invalid\n", getpid());
    sem_post(&barrier->var_mutex);
    //STAMP("%d, barrier invalid\n");
    errno = EBADF;
    return res;
  }
  barrier->waiters_cnt++;
  ////STAMP("waiters_cnt %d, val %d\n", barrier->waiters_cnt, barrier->val);
  if(barrier->waiters_cnt == barrier->val){
    //STAMP("%d, post var_mutex unlock\n", getpid());
    barrier->waiters_cnt--;
    sem_post(&barrier->var_mutex);
    for(unsigned int i = 1; i < barrier->val; ++i){
      sem_post(&barrier->b_sem); //(val - 1) post
    }
  }
  else{//release var mutex and wait on barrier semaphore
    //STAMP("%d, post var_mutex, wait b_sem\n", getpid());
    sem_post(&barrier->var_mutex);
    sem_wait(&barrier->b_sem);
    //STAMP("%d, b_sem unlocked\n", getpid());
    if(!barrier->valid){
      //STAMP("%d, Barrier invalidated\n", getpid());
      errno = ECANCELED;
      return res;
    }
  }
  sem_wait(&barrier->var_mutex);
  barrier->waiters_cnt--;
  sem_post(&barrier->var_mutex);

  //STAMP("%d, return res\n", getpid());
  return res;
}


void carma_ipcs::free_barrier(unsigned int id){
  std::map<unsigned int, sh_barrier *>::iterator it;
  it = barriers.find(id);
  if(it == barriers.end()){ //not found in map
    return;
  }
  sh_barrier *barrier = barriers[id];
  barriers.erase(it);
  //STAMP("%d, wait for var_mutex\n", getpid());
  sem_wait(&barrier->var_mutex);
  barrier->nb_proc--;
  barrier->valid = false;
  if(barrier->waiters_cnt){
    for(unsigned int i = 0; i < barrier->waiters_cnt; ++i)
      sem_post(&barrier->b_sem);
  }
  if(barrier->nb_proc == 0){
    sem_destroy(&barrier->var_mutex);
    sem_destroy(&barrier->b_sem);
    free_shm(barrier->name, barrier, sizeof(sh_barrier));
  }
  else{
    sem_post(&barrier->var_mutex);
    close_shm(barrier, sizeof(sh_barrier));
  }
}
