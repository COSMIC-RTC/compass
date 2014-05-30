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

/**
 * \brief a simple timer for CUDA kernel.
 */
class carma_timer
{
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

  void start() {
    cudaEventRecord(startEv, 0);
  }

  void reset() {
    total_time = 0.0;
  }

  /** stop timer and accumulate time in seconds */
  void stop() {
    float gpuTime;
    cudaEventRecord(stopEv, 0);
    cudaEventSynchronize(stopEv);
    cudaEventElapsedTime(&gpuTime, startEv, stopEv);
    total_time += (double)1e-3*gpuTime;
  }
    
  /** return elapsed time in seconds (as record in total_time) */
  double elapsed() {
    return total_time;
  }

};
  


#endif // CARMA_TIMER_H_
