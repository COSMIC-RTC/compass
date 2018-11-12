/**
 * \file carma_timer.h
 * \brief A simple timer class for CUDA based on events.
 *
 * \author Pierre Kestener
 * \date 30 Oct 2010
 *
 * $Id: carma_timer.h 1783 2012-02-21 10:20:07Z pkestene $
 */
#ifndef CARMA_TIMER_H_
#define CARMA_TIMER_H_

#include <carma_context.h>
#include <carma_obj.h>

/**
 * \brief a simple timer for CUDA kernel.
 */
class carma_timer {
 protected:
  cudaEvent_t startEv, stopEv;
  double total_time;

 public:
  carma_timer() {
    cudaEventCreate(&startEv);
    cudaEventCreate(&stopEv);
    total_time = 0.0;
  }

  ~carma_timer() {
    cudaEventDestroy(startEv);
    cudaEventDestroy(stopEv);
  }

  void start() { cudaEventRecord(startEv, 0); }

  void reset() { total_time = 0.0; }

  /** stop timer and accumulate time in seconds */
  void stop() {
    float gpuTime;
    cudaEventRecord(stopEv, 0);
    cudaEventSynchronize(stopEv);
    cudaEventElapsedTime(&gpuTime, startEv, stopEv);
    total_time += (double)1e-3 * gpuTime;
  }

  /** return elapsed time in seconds (as record in total_time) */
  double elapsed() { return total_time; }
};

void getClockCount(long long int *clockCounter);
void getClockCount(double *, long long int *clockCounter, double GPUfreq);
class carma_clock {
 public:
  carma_obj<double> *timeBuffer;
  double GPUfreq;
  long cc;
  long long int *clockCounter;

  carma_clock(carma_context *context, int i) {
    cudaDeviceProp cdp;
    cudaGetDeviceProperties(&cdp, context->get_activeDevice());
    GPUfreq = cdp.clockRate * 1000;
    long dims[2] = {1, i};
    timeBuffer = new carma_obj<double>(context, dims);
    cc = 0;
    cudaMalloc(&clockCounter, sizeof(long long int));
  }

  ~carma_clock() {
    delete timeBuffer;
    cudaFree(clockCounter);
  }

  void tic() { getClockCount(clockCounter); }

  void toc() {
    getClockCount(timeBuffer->getDataAt(cc), clockCounter, GPUfreq);
    cc++;
    if (cc >= timeBuffer->getNbElem()) cc = 0;
  }
};
#endif  // CARMA_TIMER_H_
