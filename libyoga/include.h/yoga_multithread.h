#ifndef YOGA_MULTITHREAD_H
#define YOGA_MULTITHREAD_H


//Simple portable thread library.

//POSIX threads.
#include <pthread.h>

typedef pthread_t yoga_thread;
typedef void *(*YOGAT_routine)(void *);

#define YOGAT_THREADPROC void*
#define YOGAT_THREADEND return 0

struct yoga_thread_barrier
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
  yoga_thread yoga_start_thread(YOGAT_routine func, void *data);
//Wait for thread to finish
  void yoga_end_thread(yoga_thread thread);
//Destroy thread
  void yoga_destroy_thread(yoga_thread thread);
//Wait for multiple threads
  void yoga_wait4thread(const yoga_thread *threads, int num);
//Create barrier.
  yoga_thread_barrier yoga_create_barrier(int releaseCount);
//Increment barrier. (excution continues)
  void yoga_increment_barrier(yoga_thread_barrier *barrier);
//Wait for barrier release.
  void yoga_wait4barrier(yoga_thread_barrier *barrier);
//Destory barrier
  void yoga_destroy_barrier(yoga_thread_barrier *barrier);

#ifdef __cplusplus
} //extern "C"
#endif

#endif //YOGA_MULTITHREAD_H
