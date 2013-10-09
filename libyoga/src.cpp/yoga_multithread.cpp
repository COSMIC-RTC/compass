#include <yoga_multithread.h>
#include <iostream>

using namespace std;

yoga_thread yoga_start_thread(YOGAT_routine func, void *data)
{
    pthread_t thread;
    pthread_create(&thread, NULL, func, data);
    return thread;
}

//Wait for thread to finish
void yoga_end_thread(yoga_thread thread)
{
    pthread_join(thread, NULL);
}

//Destroy thread
void yoga_destroy_thread(yoga_thread thread)
{
    pthread_cancel(thread);
}

//Wait for multiple threads
void yoga_wait4thread(const yoga_thread *threads, int num)
{
    for (int i = 0; i < num; i++) {
        yoga_end_thread(threads[i]);
    }
}

//Create barrier.
yoga_thread_barrier yoga_create_barrier(int releaseCount)
{
    yoga_thread_barrier barrier;

    barrier.count = 0;
    barrier.releaseCount = releaseCount;

    pthread_mutex_init(&barrier.mutex, 0);
    pthread_cond_init(&barrier.conditionVariable,0);

    return barrier;
}

//Increment barrier. (excution continues)
void yoga_increment_barrier(yoga_thread_barrier *barrier)
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
void yoga_wait4barrier(yoga_thread_barrier *barrier)
{
    pthread_mutex_lock(&barrier->mutex);

    while (barrier->count < barrier->releaseCount)
    {
        pthread_cond_wait(&barrier->conditionVariable, &barrier->mutex);
    }

    pthread_mutex_unlock(&barrier->mutex);
}

//Destory barrier
void yoga_destroy_barrier(yoga_thread_barrier *barrier)
{
    pthread_mutex_destroy(&barrier->mutex);
    pthread_cond_destroy(&barrier->conditionVariable);
}

