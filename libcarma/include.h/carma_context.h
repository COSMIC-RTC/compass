// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      carma_context.h
//! \ingroup   libcarma
//! \class     CarmaContext
//! \brief     this class provides the context in which CarmaObj are created
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.4
//! \date      2022/01/24

#ifndef _CARMA_CONTEXT_H_
#define _CARMA_CONTEXT_H_

#include <cuda_runtime_api.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector_types.h>
#include <memory>
#include <string>
#include <vector>

#include <cublas_v2.h>
#include <cusparse_v2.h>
#include <cusolverDn.h>
#include "carma_utils.h"

class CarmaDevice {
 protected:
  int id;
  cudaDeviceProp properties;
  float compute_perf;
  float cores_per_sm;
  size_t free_mem;
  size_t total_mem;

  cublasHandle_t cublas_handle;
  cusparseHandle_t cusparse_handle;
  cusolverDnHandle_t cusolver_handle;
  cudaStream_t main_stream;

 public:
  CarmaDevice(int devid);
  int set_cublas_math_mode(bool tensor);
  // CarmaDevice(const CarmaDevice& device);
  ~CarmaDevice();
  cudaStream_t get_stream() { return main_stream; }
  int get_id() { return id; }
  cudaDeviceProp get_properties() { return properties; }
  float get_compute_perf() { return compute_perf; }
  float get_cores_per_sm() { return cores_per_sm; }
  bool is_gpu_capable_p2p() { return (bool)(properties.major >= 2); }

  std::string get_name() { return properties.name; }

  size_t get_total_mem() { return total_mem; }
  size_t get_free_mem() {
    carma_safe_call(cudaMemGetInfo(&total_mem, &total_mem));
    return total_mem;
  }

  cublasHandle_t get_cublas_handle() { return cublas_handle; }
  cusparseHandle_t get_cusparse_handle() { return cusparse_handle; }
  cusolverDnHandle_t get_cusolver_handle() { return cusolver_handle; }
};

#define set_active_device(new_device, silent) \
  _set_active_device(new_device, silent, __FILE__, __LINE__)
#define set_active_device_force(new_device, silent) \
  _set_active_device_force(new_device, silent, __FILE__, __LINE__)
#define set_active_deviceForCpy(new_device, silent) \
  _set_active_device_for_copy(new_device, silent, __FILE__, __LINE__)

class CarmaContext {
 private:
  int ndevice;
  std::vector<CarmaDevice*> devices;
  int active_device;
  int** can_access_peer;

  /// singleton context
  static std::shared_ptr<CarmaContext> s_instance;

  CarmaContext();
  CarmaContext(int num_device);
  CarmaContext(int nb_devices, int32_t* devices);

  CarmaContext& operator=(const CarmaContext&) { return *s_instance; }
  CarmaContext(const CarmaContext&)
      : ndevice(0), active_device(0), can_access_peer(nullptr) {}

  void init_context(const int nb_devices, int32_t* devices_id);

 public:
  ~CarmaContext();

  static CarmaContext& instance_1gpu(int num_device);
  static CarmaContext& instance_ngpu(int nb_devices, int32_t* devices_id);
  static CarmaContext& instance();

  int get_ndevice() { return ndevice; }
  CarmaDevice* get_device(int dev) { return devices[dev]; }
  int get_active_device() { return active_device; }
  int get_active_real_device() { return devices[active_device]->get_id(); }

  int get_cuda_runtime_get_version() {
    int runtime_version;
    carma_safe_call(cudaRuntimeGetVersion(&runtime_version));
    return runtime_version;
  }

  int get_cuda_driver_get_version() {
    int driver_version;
    carma_safe_call(cudaRuntimeGetVersion(&driver_version));
    return driver_version;
  }

  std::string get_device_name(int device);
  std::string get_device_info(int device);
  std::string get_device_mem_info(int device);

  inline int _set_active_device_for_copy(int new_device, int silent,
                                     std::string file, int line) {
    if (new_device > ndevice) return -1;
    return (can_access_peer[active_device][new_device] != 1)
               ? _set_active_device(new_device, silent, file, line)
               : active_device;
  }
  inline int _set_active_device(int new_device, int silent, std::string file,
                               int line) {
    return (this->active_device != new_device)
               ? _set_active_device_force(new_device, silent, file, line)
               : active_device;
  }
  int _set_active_device_force(int new_device, int silent, std::string file,
                             int line);
  int get_max_gflops_device_id();
  cublasHandle_t get_cublas_handle() { return get_cublas_handle(active_device); }
  cusparseHandle_t get_cusparse_handle() {
    return get_cusparse_handle(active_device);
  }

  cublasHandle_t get_cublas_handle(int device) {
    return devices[device]->get_cublas_handle();
  }
  cusparseHandle_t get_cusparse_handle(int device) {
    return devices[device]->get_cusparse_handle();
  }
  bool can_p2p(int dev1, int dev2) { return can_access_peer[dev1][dev2]; }

  std::string magma_info();
};

/// from /usr/local/cuda/samples/common/inc/helper_cuda.h
inline int convert_sm_version2cores(int major, int minor) {
  // Defines for GPU Architecture types (using the SM version to determine the #
  // of cores per SM
  typedef struct {
    int SM;  // 0xMm (hexidecimal notation), M = SM Major version, and m = SM
             // minor version
    int Cores;
  } sSMtoCores;

  sSMtoCores nGpuArchCoresPerSM[] = {
      {0x20, 32},   // Fermi Generation (SM 2.0) GF100 class
      {0x21, 48},   // Fermi Generation (SM 2.1) GF10x class
      {0x30, 192},  // Kepler Generation (SM 3.0) GK10x class
      {0x32, 192},  // Kepler Generation (SM 3.2) GK10x class
      {0x35, 192},  // Kepler Generation (SM 3.5) GK11x class
      {0x37, 192},  // Kepler Generation (SM 3.7) GK21x class
      {0x50, 128},  // Maxwell Generation (SM 5.0) GM10x class
      {0x52, 128},  // Maxwell Generation (SM 5.2) GM20x class
      {0x53, 128},  // Maxwell Generation (SM 5.3) GM20x class
      {0x60, 64},   // Pascal Generation (SM 6.0) GP100 class
      {0x61, 128},  // Pascal Generation (SM 6.1) GP10x class
      {0x62, 128},  // Pascal Generation (SM 6.2) GP10x class
      {0x70, 64},   // Volta Generation (SM 7.0) GV100 class
      {0x72, 64},   // Volta Generation (SM 7.2) GV10B class
      {0x75, 64},   // Turing Generation (SM 7.5) TU100 class
      {0x80, 64},   // Ampere Generation (SM 8.0) GA102 class
      {0x86, 64},   // Ampere Generation (SM 8.6) GA104 class
      {-1, -1}};

  int index = 0;

  while (nGpuArchCoresPerSM[index].SM != -1) {
    if (nGpuArchCoresPerSM[index].SM == ((major << 4) + minor)) {
      return nGpuArchCoresPerSM[index].Cores;
    }

    index++;
  }

  // If we don't find the values, we default use the previous one to run
  // properly
  printf(
      "MapSMtoCores for SM %d.%d is undefined.  Default to use %d Cores/SM\n",
      major, minor, nGpuArchCoresPerSM[index - 1].Cores);
  return nGpuArchCoresPerSM[index - 1].Cores;
}
#endif  // _CARMA_CONTEXT_H_
