// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      carma_timer.hpp
//! \ingroup   libcarma
//! \class     CarmaTimer
//! \brief     this class provides the timer features to CarmaObj
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24


#ifndef CARMA_TIMER_H_
#define CARMA_TIMER_H_

#include <carma_context.hpp>
#include <carma_obj.hpp>

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

void get_clock_count(int64_t *clock_counter);
void get_clock_count(double *, int64_t *clock_counter, double gpu_freq);
class CarmaClock {
 public:
  CarmaObj<double> *time_buffer;
  double gpu_freq;
  int64_t cc;
  int64_t *clock_counter;

  CarmaClock(CarmaContext *context, int32_t i) {
    cudaDeviceProp cdp;
    cudaGetDeviceProperties(&cdp, context->get_active_device());
    gpu_freq = cdp.clockRate * 1000;
    int64_t dims[2] = {1, i};
    time_buffer = new CarmaObj<double>(context, dims);
    cc = 0;
    cudaMalloc(&clock_counter, sizeof(int64_t));
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
