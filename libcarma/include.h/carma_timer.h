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

//! \file      carma_timer.h
//! \ingroup   libcarma
//! \class     CarmaTimer
//! \brief     this class provides the timer features to CarmaObj
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.1.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License


#ifndef CARMA_TIMER_H_
#define CARMA_TIMER_H_

#include <carma_context.h>
#include <carma_obj.h>

/**
 * \brief a simple timer for CUDA kernel.
 */
class CarmaTimer {
 protected:
  cudaEvent_t start_event, stop_event;
  cudaStream_t stream = 0;
  double total_time;

 public:
  CarmaTimer() {
    carma_safe_call(cudaEventCreate(&start_event));
    carma_safe_call(cudaEventCreate(&stop_event));
    total_time = 0.0;
  }

  ~CarmaTimer() {
    cudaEventDestroy(start_event);
    cudaEventDestroy(stop_event);
  }

  void start() { carma_safe_call(cudaEventRecord(start_event, stream)); }

  void reset() { total_time = 0.0; }

  /** stop timer and accumulate time in seconds */
  void stop() {
    float gpuTime;
    carma_safe_call(cudaEventRecord(stop_event, stream));
    carma_safe_call(cudaEventSynchronize(stop_event));
    carma_safe_call(cudaEventElapsedTime(&gpuTime, start_event, stop_event));
    total_time += (double)1e-3 * gpuTime;
  }

  void set_stream(cudaStream_t newStream) { stream = newStream; }
  /** return elapsed time in seconds (as record in total_time) */
  double elapsed() { return total_time; }
};

void get_clock_count(long long int *clock_counter);
void get_clock_count(double *, long long int *clock_counter, double gpu_freq);
class CarmaClock {
 public:
  CarmaObj<double> *time_buffer;
  double gpu_freq;
  long cc;
  long long int *clock_counter;

  CarmaClock(CarmaContext *context, int i) {
    cudaDeviceProp cdp;
    cudaGetDeviceProperties(&cdp, context->get_active_device());
    gpu_freq = cdp.clockRate * 1000;
    long dims[2] = {1, i};
    time_buffer = new CarmaObj<double>(context, dims);
    cc = 0;
    cudaMalloc(&clock_counter, sizeof(long long int));
  }

  ~CarmaClock() {
    delete time_buffer;
    cudaFree(clock_counter);
  }

  void tic() { get_clock_count(clock_counter); }

  void toc() {
    get_clock_count(time_buffer->get_data_at(cc), clock_counter, gpu_freq);
    cc++;
    if (cc >= time_buffer->get_nb_elements()) cc = 0;
  }
};
#endif  // CARMA_TIMER_H_
