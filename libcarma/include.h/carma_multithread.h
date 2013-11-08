#ifndef CARMA_MULTITHREAD_H
#define CARMA_MULTITHREAD_H


//Simple portable thread library.

//POSIX threads.
#include <pthread.h>

typedef pthread_t carma_thread;
typedef void *(*CARMAT_routine)(void *);

#define CARMAT_THREADPROC void*
#define CARMAT_THREADEND return 0

struct carma_thread_barrier
{
    pthread_mutex_t mutex;
    pthread_cond_t conditionVariable;
    int releaseCount;
    int count;
};

#ifdef __cplusplus
extern "C" {
#endif

//Create thread
  carma_thread carma_start_thread(CARMAT_routine func, void *data);
//Wait for thread to finish
  void carma_end_thread(carma_thread thread);
//Destroy thread
  void carma_destroy_thread(carma_thread thread);
//Wait for multiple threads
  void carma_wait4thread(const carma_thread *threads, int num);
//Create barrier.
  carma_thread_barrier carma_create_barrier(int releaseCount);
//Increment barrier. (excution continues)
  void carma_increment_barrier(carma_thread_barrier *barrier);
//Wait for barrier release.
  void carma_wait4barrier(carma_thread_barrier *barrier);
//Destory barrier
  void carma_destroy_barrier(carma_thread_barrier *barrier);

#ifdef __cplusplus
} //extern "C"
#endif

#endif //CARMA_MULTITHREAD_H
