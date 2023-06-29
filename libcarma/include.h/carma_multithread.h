// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      carma_multithread.h
//! \ingroup   libcarma
//! \brief     this fle provides the multithread features to CarmaObj
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.4
//! \date      2022/01/24

#ifndef CARMA_MULTITHREAD_H
#define CARMA_MULTITHREAD_H

// Simple portable thread library.

// POSIX threads.
#include <pthread.h>

typedef pthread_t carma_thread;
typedef void *(*carma_routine)(void *);

#define CARMAT_THREADPROC void *
#define CARMAT_THREADEND return 0

struct CarmaThreadBarrier {
  pthread_mutex_t mutex;
  pthread_cond_t conditionVariable;
  int releaseCount;
  int count;
};

#ifdef __cplusplus
extern "C" {
#endif

// Create thread
carma_thread carma_start_thread(carma_routine func, void *data);
// Wait for thread to finish
void carma_end_thread(carma_thread thread);
// Destroy thread
void carma_destroy_thread(carma_thread thread);
// Wait for multiple threads
void carma_wait4thread(const carma_thread *threads, int num);
// Create barrier.
CarmaThreadBarrier carma_create_barrier(int releaseCount);
// Increment barrier. (excution continues)
void carma_increment_barrier(CarmaThreadBarrier *barrier);
// Wait for barrier release.
void carma_wait4barrier(CarmaThreadBarrier *barrier);
// Destory barrier
void carma_destroy_barrier(CarmaThreadBarrier *barrier);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // CARMA_MULTITHREAD_H
