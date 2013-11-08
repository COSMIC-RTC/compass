#include <carma_multithread.h>
#include <iostream>

using namespace std;

carma_thread carma_start_thread(CARMAT_routine func, void *data)
{
    pthread_t thread;
    pthread_create(&thread, NULL, func, data);
    return thread;
}

//Wait for thread to finish
void carma_end_thread(carma_thread thread)
{
    pthread_join(thread, NULL);
}

//Destroy thread
void carma_destroy_thread(carma_thread thread)
{
    pthread_cancel(thread);
}

//Wait for multiple threads
void carma_wait4thread(const carma_thread *threads, int num)
{
    for (int i = 0; i < num; i++) {
        carma_end_thread(threads[i]);
    }
}

//Create barrier.
carma_thread_barrier carma_create_barrier(int releaseCount)
{
    carma_thread_barrier barrier;

    barrier.count = 0;
    barrier.releaseCount = releaseCount;

    pthread_mutex_init(&barrier.mutex, 0);
    pthread_cond_init(&barrier.conditionVariable,0);

    return barrier;
}

//Increment barrier. (excution continues)
void carma_increment_barrier(carma_thread_barrier *barrier)
{
    int myBarrierCount;
    pthread_mutex_lock(&barrier->mutex);
    myBarrierCount = ++barrier->count;
    pthread_mutex_unlock(&barrier->mutex);

    if (myBarrierCount >=barrier->releaseCount)
    {
        pthread_cond_signal(&barrier->conditionVariable);
    }
}

//Wait for barrier release.
void carma_wait4barrier(carma_thread_barrier *barrier)
{
    pthread_mutex_lock(&barrier->mutex);

    while (barrier->count < barrier->releaseCount)
    {
        pthread_cond_wait(&barrier->conditionVariable, &barrier->mutex);
    }

    pthread_mutex_unlock(&barrier->mutex);
}

//Destory barrier
void carma_destroy_barrier(carma_thread_barrier *barrier)
{
    pthread_mutex_destroy(&barrier->mutex);
    pthread_cond_destroy(&barrier->conditionVariable);
}

