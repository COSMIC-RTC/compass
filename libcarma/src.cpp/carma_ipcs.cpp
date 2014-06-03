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
  //else, everything's OK
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
  sprintf(name, "/cipcs_cubuffer_%d_in", id);
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
  buffer->nb_proc = 1;
  strcpy(buffer->name, name);
  buffer->isBoard = isBoard;

  sem_post(&buffer->mutex); //mutex = 1;

  buffers[id] = buffer;

  errno = 0;
  res = EXIT_SUCCESS;
  return res;
}


int carma_ipcs::map_transfer_shm(unsigned int id){
  int res = EXIT_FAILURE;
  std::map<unsigned int, sh_buffer *>::iterator it;
    it = buffers.find(id);
    if(it == buffers.end()){ //not found in map
      sh_buffer *buffer;
      char name[NAME_MAX+1];
      sprintf(name, "/cipcs_cubuffer_%d", id);
      if((buffer = (sh_buffer *)get_shm(name)) == (void *) -1){ //not found in shm
        errno = EINVAL;
        return res;
      }
      else{ //found in shm
        if(sem_wait(&buffer->var_mutex))//probably destroyed while waiting on it, i.e. no more sharing
          return res;
        if(buffer->owner == getpid()){ //owner try to get handler, illegal operation
          sem_post(&buffer->var_mutex);
          errno = EPERM;
          return res;
        }
        //else, add it to the local map
        buffer->nb_proc++;
        sem_post(&buffer->var_mutex);
        buffers[id] = buffer;
      }
    }
}


int carma_ipcs::write_gpu(CUdeviceptr src, void *dst, size_t bsize){
  int res = EXIT_FAILURE;
  return res;
}


int carma_ipcs::read_gpu(void *src, CUdeviceptr dst, size_t bsize){
  int res = EXIT_FAILURE;
  return res;
  //todo
}


int carma_ipcs::write_transfer_shm(unsigned int id, void *src, size_t bsize, bool gpuBuffer){
  int res = EXIT_FAILURE;
  return res;
  //todo
}


int carma_ipcs::read_transfer_shm(unsigned int id, void *dst, size_t bsize, bool gpuBuffer){
  int res = EXIT_FAILURE;
  return res;
  //todo
}

int unmap_tranfer_shm(unsigned int id){

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
  if(!buffer->nb_proc){
    if(buffer->isBoard)
      sem_destroy(&buffer->wait_pub);
    sem_destroy(&buffer->mutex);
    char in_name[NAME_MAX+1];
    sprintf(in_name, "%s_in", buffer->name);

    free_shm(in_name, buffer->p_shm, buffer->size);
    free_shm(buffer->name, buffer, sizeof(sh_buffer));
  }
  else{
    sem_post(&buffer->var_mutex);
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
  barrier->proc_cnt = 0;
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
  it = barriers.find(id);
  if(it == barriers.end()){ //not found in map
    //STAMP("barrier not found in map\n");
    sh_barrier * barrier;
    char name[NAME_MAX+1];

    sprintf(name, "/cipcs_cubarrier_%d", id);
    if((barrier = (sh_barrier *)get_shm(name)) == (void *) -1){ //not found in shm
      //STAMP("barrier not found in shm\n");
      errno = EINVAL;
      return res;
    }
    //else, found in shm, add to map
    barriers[id] = barrier;
    //STAMP("barrier added to map\n");

  }

  //access barrier
  sh_barrier *barrier = barriers[id];
  //STAMP("wait for var_mutex\n");
  if(!barrier->valid && sem_wait(&barrier->var_mutex) == -1) //probably invalided, then destroyed, probability low...
    return res;
  if(!barrier->valid){
    //STAMP("barrier invalid\n");
    free_barrier(id);
    errno = EAGAIN;
    return res;
  }
  barrier->proc_cnt++;
  if(barrier->proc_cnt == barrier->val){
    for(unsigned int i = 1; i < barrier->val; ++i)
      sem_post(&barrier->b_sem); //(val - 1) post
    barrier->proc_cnt = 0; //reinit barrier
    sem_post(&barrier->var_mutex);
  }
  else{//release var mutex and wait on barrier semaphore
    sem_post(&barrier->var_mutex);
    sem_wait(&barrier->b_sem);
    if(!barrier->valid){
      errno = ECANCELED;
      return res;
    }
  }

  errno = 0;
  res = EXIT_SUCCESS;
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
//  STAMP("wait for var_mutex\n");
  sem_wait(&barrier->var_mutex);
  barrier->valid = false;
  if(!barrier->proc_cnt){
    sem_destroy(&barrier->var_mutex);
    sem_destroy(&barrier->b_sem);
    free_shm(barrier->name, barrier, sizeof(sh_barrier));
  }
  else{
    for(unsigned int i = 0; i < barrier->proc_cnt; ++i)
      sem_post(&barrier->b_sem);
    barrier->proc_cnt = 0;
    sem_post(&barrier->var_mutex);
    close_shm(barrier, sizeof(sh_barrier));
  }
}
