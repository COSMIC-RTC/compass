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

//! \file      carma_multithread.h
//! \ingroup   libcarma
//! \brief     this fle provides the multithread features to CarmaObj
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.2.1
//! \date      2022/01/24
//! \copyright GNU Lesser General Public License

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
