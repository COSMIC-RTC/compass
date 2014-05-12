/*
 * prana_rtc.h
 *
 *  Created on: May 9, 2014
 *      Author: sevin
 */

#ifndef PRANA_RTC_H_
#define PRANA_RTC_H_

#include<carma_multithread.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <vector>
#include <vector_types.h>

using std::vector;

class prana_rtc {
public:
  prana_rtc(int nGPUs, CUcontext *ctx,
      cublasHandle_t cublas_handle, int nPixImgX, int nPixImgY, int nSlp,
      int nCmd, int nSppX, int nSppY, short *h_iSppValid, float *h_matcom);
  virtual ~prana_rtc();

  void set_affinity(int cpu);

  void start_rtc();
  void stop_rtc();

  void set_matcom(CUdeviceptr d_matcom, CUcontext ctx);
  void set_image(CUdeviceptr d_image, CUcontext ctx);
  void get_commands(CUdeviceptr d_commands, CUcontext ctx);

private:

  int nGPUs;
  CUcontext *ctx;
  int ctx_use;

  cublasStatus_t status_cublas;
  cublasHandle_t cublas_handle;

  CUresult err;
  CUmodule cuModule_centroider;
  CUmodule cuModule_acquisim;
  CUmodule cuModule_sync;

  CUfunction reduce2, centroidx, centroidy, bcube_krnl;

  int nPixImgX;
  int nPixImgY;
  int nSlp;
  int nCmd;
  int nSppX;
  int nSppY;

  CUdeviceptr d_image_full;
  CUdeviceptr d_image_cube;
  CUdeviceptr d_iSppValid;
  CUdeviceptr d_com_full;

  vector<CUdeviceptr> d_image;
  vector<CUdeviceptr> d_cmat;
  vector<CUdeviceptr> d_slopes;
  vector<CUdeviceptr> d_subsum;
  vector<CUdeviceptr> d_com;

  CUevent start, stop;

  carma_thread_barrier start_compute, stop_compute;
  carma_thread rtc_loop;

  bool running;
  void *run();
  static void *run_helper(void * data);
};

#endif /* PRANA_RTC_H_ */
