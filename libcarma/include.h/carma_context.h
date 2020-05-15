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

//! \file      carma_context.h
//! \ingroup   libcarma
//! \class     carma_context
//! \brief     this class provides the context in which carma_obj are created
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.4.1
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#ifndef _CARMA_CONTEXT_H_
#define _CARMA_CONTEXT_H_

#include <cuda_runtime_api.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector_types.h>
#include <memory>
#include <string>
#include <vector>

#include "carma_cublas.h"
#include "carma_cusparse.h"
#include "carma_utils.h"

class carma_device {
 protected:
  int id;
  cudaDeviceProp properties;
  float compute_perf;
  float cores_per_sm;
  bool p2p_activate;
  size_t freeMem;
  size_t totalMem;

  cublasHandle_t cublasHandle;
  cusparseHandle_t cusparseHandle;
  cudaStream_t mainStream;

 public:
  carma_device(int devid);
  int set_cublas_math_mode(bool tensor);
  // carma_device(const carma_device& device);
  ~carma_device();
  cudaStream_t get_stream() { return mainStream; }
  int get_id() { return id; }
  cudaDeviceProp get_properties() { return properties; }
  float get_compute_perf() { return compute_perf; }
  float get_cores_per_sm() { return cores_per_sm; }
  bool isGPUCapableP2P() { return (bool)(properties.major >= 2); }

  bool isP2P_active() { return p2p_activate; }

  std::string getName() { return properties.name; }

  size_t getTotalMem() { return totalMem; }
  size_t getFreeMem() {
    carmaSafeCall(cudaMemGetInfo(&freeMem, &totalMem));
    return freeMem;
  }

  cublasHandle_t get_cublasHandle() { return cublasHandle; }
  cusparseHandle_t get_cusparseHandle() { return cusparseHandle; }
};

#define set_active_device(newDevice, silent) \
  _set_active_device(newDevice, silent, __FILE__, __LINE__)
#define set_active_device_force(newDevice, silent) \
  _set_active_device_force(newDevice, silent, __FILE__, __LINE__)
#define set_active_deviceForCpy(newDevice, silent) \
  _set_active_deviceForCpy(newDevice, silent, __FILE__, __LINE__)

class carma_context {
 private:
  int ndevice;
  std::vector<carma_device*> devices;
  int active_device;
  int** can_access_peer;

  /// singleton context
  static std::shared_ptr<carma_context> s_instance;

  carma_context();
  carma_context(int num_device);
  carma_context(int nb_devices, int32_t* devices);

  carma_context& operator=(const carma_context&) { return *s_instance; }
  carma_context(const carma_context&)
      : ndevice(0), active_device(0), can_access_peer(nullptr) {}

  void init_context(const int nb_devices, int32_t* devices_id);

 public:
  ~carma_context();

  static carma_context& instance_1gpu(int num_device);
  static carma_context& instance_ngpu(int nb_devices, int32_t* devices_id);
  static carma_context& instance();

  int get_ndevice() { return ndevice; }
  carma_device* get_device(int dev) { return devices[dev]; }
  int get_active_device() { return active_device; }
  int get_activeRealDevice() { return devices[active_device]->get_id(); }

  int get_cudaRuntimeGetVersion() {
    int runtimeVersion;
    carmaSafeCall(cudaRuntimeGetVersion(&runtimeVersion));
    return runtimeVersion;
  }

  int get_cudaDriverGetVersion() {
    int driverVersion;
    carmaSafeCall(cudaRuntimeGetVersion(&driverVersion));
    return driverVersion;
  }

  std::string get_DeviceName(int device);
  std::string get_DeviceInfo(int device);
  std::string get_DeviceMemInfo(int device);

  inline int _set_active_deviceForCpy(int newDevice, int silent,
                                     std::string file, int line) {
    if (newDevice > ndevice) return -1;
    return (can_access_peer[active_device][newDevice] != 1)
               ? _set_active_device(newDevice, silent, file, line)
               : active_device;
  }
  inline int _set_active_device(int newDevice, int silent, std::string file,
                               int line) {
    return (this->active_device != newDevice)
               ? _set_active_device_force(newDevice, silent, file, line)
               : active_device;
  }
  int _set_active_device_force(int newDevice, int silent, std::string file,
                             int line);
  int get_maxGflopsDeviceId();
  cublasHandle_t get_cublasHandle() { return get_cublasHandle(active_device); }
  cusparseHandle_t get_cusparseHandle() {
    return get_cusparseHandle(active_device);
  }

  cublasHandle_t get_cublasHandle(int device) {
    return devices[device]->get_cublasHandle();
  }
  cusparseHandle_t get_cusparseHandle(int device) {
    return devices[device]->get_cusparseHandle();
  }
  bool canP2P(int dev1, int dev2) { return can_access_peer[dev1][dev2]; }

  std::string magma_info();
};

/// from /usr/local/cuda/samples/common/inc/helper_cuda.h
inline int ConvertSMVer2Cores(int major, int minor) {
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
