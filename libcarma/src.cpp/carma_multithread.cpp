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

//! \file      carma_multithread.cpp
//! \ingroup   libcarma
//! \brief     this fle provides the multithread features to CarmaObj
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.2.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

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
void carma_wait4thread(const carma_thread *threads, int num) {
  for (int i = 0; i < num; i++) {
    carma_end_thread(threads[i]);
  }
}

// Create barrier.
CarmaThreadBarrier carma_create_barrier(int releaseCount) {
  CarmaThreadBarrier barrier;

  barrier.count = 0;
  barrier.releaseCount = releaseCount;

  pthread_mutex_init(&barrier.mutex, 0);
  pthread_cond_init(&barrier.conditionVariable, 0);

  return barrier;
}

// Increment barrier. (excution continues)
void carma_increment_barrier(CarmaThreadBarrier *barrier) {
  int myBarrierCount;
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
