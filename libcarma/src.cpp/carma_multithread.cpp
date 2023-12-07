// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      carma_multithread.cpp
//! \ingroup   libcarma
//! \brief     this fle provides the multithread features to CarmaObj
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#include <carma_multithread.h>
#include <iostream>

carma_thread carma_start_thread(carma_routine func, void *data) {
  pthread_t thread;
  pthread_create(&thread, NULL, func, data);
  return thread;
}

// Wait for thread to finish
void carma_end_thread(carma_thread thread) { pthread_join(thread, NULL); }

// Destroy thread
void carma_destroy_thread(carma_thread thread) { pthread_cancel(thread); }

// Wait for multiple threads
void carma_wait4thread(const carma_thread *threads, int32_t num) {
  for (int32_t i = 0; i < num; i++) {
    carma_end_thread(threads[i]);
  }
}

// Create barrier.
CarmaThreadBarrier carma_create_barrier(int32_t releaseCount) {
  CarmaThreadBarrier barrier;

  barrier.count = 0;
  barrier.releaseCount = releaseCount;

  pthread_mutex_init(&barrier.mutex, 0);
  pthread_cond_init(&barrier.conditionVariable, 0);

  return barrier;
}

// Increment barrier. (excution continues)
void carma_increment_barrier(CarmaThreadBarrier *barrier) {
  int32_t myBarrierCount;
  pthread_mutex_lock(&barrier->mutex);
  myBarrierCount = ++barrier->count;
  pthread_mutex_unlock(&barrier->mutex);

  if (myBarrierCount >= barrier->releaseCount) {
    pthread_cond_signal(&barrier->conditionVariable);
  }
}

// Wait for barrier release.
void carma_wait4barrier(CarmaThreadBarrier *barrier) {
  pthread_mutex_lock(&barrier->mutex);

  while (barrier->count < barrier->releaseCount) {
    pthread_cond_wait(&barrier->conditionVariable, &barrier->mutex);
  }
  barrier->count = 0;
  pthread_mutex_unlock(&barrier->mutex);
}

// Destory barrier
void carma_destroy_barrier(CarmaThreadBarrier *barrier) {
  pthread_mutex_destroy(&barrier->mutex);
  pthread_cond_destroy(&barrier->conditionVariable);
}
