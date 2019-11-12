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

//! \file      carma_context.cpp
//! \ingroup   libcarma
//! \class     carma_context
//! \brief     this class provides the context in which carma_obj are created
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.3.2
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include <carma_context.h>

#ifdef USE_CULA
// CULA headers
#include <cula.hpp>
#endif  // USE_CULA

#ifdef USE_MAGMA
// MAGMA headers
#include "magma.h"
#include "magma_lapack.h"
#endif  // USE_MAGMA

carma_device::carma_device(int devid) {
  carmaSafeCall(cudaSetDevice(devid));

  // Instruct CUDA to yield its thread when waiting for results from the device.
  // This can increase latency when waiting for the device, but can increase the
  // performance of CPU threads performing work in parallel with the device. see
  // also: cudaDeviceScheduleAuto, cudaDeviceScheduleSpin,
  // cudaDeviceScheduleYield
  // carmaSafeCall(cudaSetDeviceFlags(cudaDeviceScheduleYield));

  this->id = devid;
  cudaGetDeviceProperties(&(this->properties), devid);
  this->cores_per_sm =
      ConvertSMVer2Cores(this->properties.major, this->properties.minor);
  this->compute_perf = this->properties.multiProcessorCount *
                       this->cores_per_sm * this->properties.clockRate;

  this->p2p_activate = false;

  carmaSafeCall(cudaMemGetInfo(&freeMem, &totalMem));

  carma_initCublas(&cublasHandle);
  carma_initCusparse(&cusparseHandle);

  // cudaStreamCreate(&mainStream);
  //  cusparsePointerMode_t mode;
  //  cusparseGetPointerMode(cusparseHandle, &mode);
  //  DEBUG_TRACE("%d\n", mode);

  // DEBUG_TRACE("done\n");
}
int carma_device::set_cublas_math_mode(bool tensor) {
  if (tensor)
    carma_checkCublasStatus(
        cublasSetMathMode(cublasHandle, CUBLAS_TENSOR_OP_MATH));
  else
    carma_checkCublasStatus(
        cublasSetMathMode(cublasHandle, CUBLAS_DEFAULT_MATH));

  return EXIT_SUCCESS;
}

carma_device::~carma_device() {
  carma_shutdownCublas(cublasHandle);
  carma_shutdownCusparse(cusparseHandle);
  // cudaStreamDestroy(mainStream);

  this->id = -1;
}

std::shared_ptr<carma_context> carma_context::s_instance;

carma_context &carma_context::instance_1gpu(int num_device) {
  if (!carma_context::s_instance) {
    carma_context::s_instance =
        std::shared_ptr<carma_context>(new carma_context(num_device));
  }
  return *carma_context::s_instance;
}

carma_context &carma_context::instance_ngpu(int nb_devices,
                                            int32_t *devices_id) {
  if (!carma_context::s_instance) {
    carma_context::s_instance = std::shared_ptr<carma_context>(
        new carma_context(nb_devices, devices_id));
  }
  return *carma_context::s_instance;
}

carma_context &carma_context::instance() {
  if (!carma_context::s_instance) {
    carma_context::s_instance =
        std::shared_ptr<carma_context>(new carma_context());
  }
  return *carma_context::s_instance;
}

carma_context::carma_context(int num_device) {
  can_access_peer = nullptr;
  this->activeDevice = -1;
  this->ndevice = -1;

  int devices[1];
  devices[0] = num_device;
  init_context(1, devices);
}

carma_context::carma_context() {
  carmaSafeCall(cudaGetDeviceCount(&(this->ndevice)));
  can_access_peer = nullptr;
  this->activeDevice = -1;

  if (this->ndevice == 0) {
    DEBUG_TRACE("carma_context() CUDA error: no devices supporting CUDA.");
    throw std::runtime_error(
        "carma_context() CUDA error: no devices supporting CUDA.");
  }

  int const size = this->ndevice;
  int32_t devices_id[size];

  for (int i = 0; i < size; ++i) devices_id[i] = i;
  init_context(this->ndevice, devices_id);
}

carma_context::carma_context(int nb_devices, int32_t *devices_id) {
  can_access_peer = nullptr;
  this->activeDevice = -1;
  this->ndevice = -1;

  init_context(nb_devices, devices_id);
}

void carma_context::init_context(const int nb_devices, int32_t *devices_id) {
  // TODO(all) : why seed is initialized here ?
  srandom(1234);
  this->activeDevice = -1;

  int n_cuda_devices = 0;

  carmaSafeCall(cudaGetDeviceCount(&n_cuda_devices));

  if (!n_cuda_devices) {
    DEBUG_TRACE("carma_context() CUDA error: no devices supporting CUDA.");
    throw std::runtime_error(
        "carma_context() CUDA error: no devices supporting CUDA.");
  }

  if (nb_devices > n_cuda_devices) {
    DEBUG_TRACE(
        "carma_context() CUDA error: not enough devices supporting CUDA. ask "
        "%d, available %d",
        nb_devices, n_cuda_devices);
    DEBUG_TRACE("carma_context() will be initialized on GPU 0 only");
    this->ndevice = 1;
  } else
    this->ndevice = nb_devices;
  int current_device = 0;

  while (current_device < this->ndevice) {
    devices.push_back(new carma_device(devices_id[current_device]));
    current_device++;
  }

  can_access_peer = new int *[this->ndevice];

  for (int i = 0; i < ndevice; i++) {
    can_access_peer[i] = new int[this->ndevice];

    for (int j = 0; j < ndevice; j++) {
      can_access_peer[i][j] = (i == j);
    }
  }

#ifdef USE_UVA

  int gpuid[this->ndevice];  // we want to find the first two GPU's that can
                             // support P2P
  int gpu_count = 0;         // GPUs that meet the criteria
  current_device = 0;

  while (current_device < this->ndevice) {
    if (devices[current_device]->isGPUCapableP2P())
      gpuid[gpu_count++] = current_device;
    current_device++;
  }

  if (gpu_count > 1) {
    bool has_uva = true;

    for (int i = 0; i < gpu_count - 1; i++) {
      has_uva &= devices[gpuid[i]]->get_properties().unifiedAddressing;

      for (int j = i + 1; j < gpu_count; j++) {
        carmaSafeCall(cudaDeviceCanAccessPeer(
            &can_access_peer[gpuid[i]][gpuid[j]], devices_id[gpuid[i]],
            devices_id[gpuid[j]]));
        carmaSafeCall(cudaDeviceCanAccessPeer(
            &can_access_peer[gpuid[j]][gpuid[i]], devices_id[gpuid[j]],
            devices_id[gpuid[i]]));

        if ((can_access_peer[gpuid[i]][gpuid[j]] == 1) &&
            (can_access_peer[gpuid[j]][gpuid[i]] == 1)) {
          printf("*** Enabling peer access between GPU%d and GPU%d... ***\n",
                 devices_id[gpuid[i]], devices_id[gpuid[j]]);
          carmaSafeCall(cudaSetDevice(devices_id[gpuid[i]]));
          carmaSafeCall(cudaDeviceEnablePeerAccess(devices_id[gpuid[j]], 0));
          carmaSafeCall(cudaSetDevice(devices_id[gpuid[j]]));
          carmaSafeCall(cudaDeviceEnablePeerAccess(devices_id[gpuid[i]], 0));
        }
      }
    }
    has_uva &=
        devices[gpuid[gpu_count - 1]]->get_properties().unifiedAddressing;

    if (has_uva) {
      printf("*** All GPUs listed can support UVA... ***\n");
    }
  }
#endif  // USE_UVA

  this->activeDevice =
      set_activeDeviceForce(0, 1);  // get_maxGflopsDeviceId(), 1);

#ifdef USE_CULA

  // CULA init
  culaStatus status = culaInitialize();

  if (status) {
    char buf[256];
    culaGetErrorInfoString(status, culaGetErrorInfo(), buf, sizeof(buf));
    printf("%s\n", buf);
  }
#endif  // USE_CULA

#ifdef USE_MAGMA

// MAGMA init
#ifdef USE_MAGMA_PATCHED
  magma_init(nb_devices, devices_id);
#else   // ifdef USE_MAGMA_PATCHED
  magma_init();
#endif  // USE_MAGMA_PATCHED

#if DEBUG

//  magma_print_environment();
#endif  // DEBUG
#endif  // USE_MAGMA

#if DEBUG
  printf("CARMA Context created @ %p\n", this);
#endif  // DEBUG
}

carma_context::~carma_context() {
  carmaSafeCall(cudaDeviceSynchronize());

#ifdef USE_CULA
  // CULA finalize
  culaShutdown();
#endif  // USE_CULA

#ifdef USE_MAGMA

  // MAGMA finalize
  magma_finalize();
#endif  // USE_MAGMA

  size_t idx = 0;

  while (this->devices.size() > 0) {
    delete this->devices.back();
    this->devices.pop_back();

    if (can_access_peer[idx] != nullptr) delete[] can_access_peer[idx];
    ++idx;
  }

  if (can_access_peer != nullptr) delete[] can_access_peer;

#if DEBUG
  printf("CARMA Context deleted @ %p\n", this);
#endif  // DEBUG
}

int carma_context::_set_activeDeviceForce(int newDevice, int silent,
                                          std::string file, int line) {
  if (newDevice < ndevice) {
    carmaSafeCall(cudaSetDevice(devices[newDevice]->get_id()));
#ifdef USE_CULA
    culaStatus status = culaSelectDevice(newDevice);

    if (status) {
      char buf[256];
      culaGetErrorInfoString(status, culaGetErrorInfo(), buf, sizeof(buf));
      printf("%s\n", buf);
    }
#endif  // USE_CULA
#if DEBUG
    silent = 0;
#endif  // DEBUG

    if (!silent) {
      std::cout << "Using device " << devices[newDevice]->get_id() << ": \""
                << devices[newDevice]->get_properties().name
                << "\" with Compute "
                << devices[newDevice]->get_properties().major << "."
                << devices[newDevice]->get_properties().minor << " capability"
                << std::endl;
    }
    activeDevice = newDevice;
  } else {
    fprintf(stderr,
            "[%s:%d] Invalid Device Id : %d, Your system has only %d CUDA "
            "capable device(s) available ",
            file.c_str(), line, newDevice, ndevice);
    std::cerr << "Leaving activeDevice to its current value : " << activeDevice
              << std::endl;
  }
  return activeDevice;
}

std::string carma_context::get_DeviceName(int device) {
  return devices[device]->get_properties().name;
}

std::string carma_context::get_DeviceInfo(int device) {
  std::stringstream buf;

  buf << "device " << device << ": \"" << devices[device]->get_properties().name
      << "\" with Compute " << devices[device]->get_properties().major << "."
      << devices[device]->get_properties().minor << " capability";
  return buf.str();
}

std::string carma_context::get_DeviceMemInfo(int device) {
  std::stringstream buf;
  size_t totalMem = devices[device]->getTotalMem() / 1024 / 1024;
  size_t usedMem = totalMem - devices[device]->getFreeMem() / 1024 / 1024;

  buf << "device " << device << ": \"" << devices[device]->get_properties().name
      << "\" memory used " << usedMem << "MB / " << totalMem << "MB ("
      << usedMem * 100. / totalMem << "%)";
  return buf.str();
}

int carma_context::get_maxGflopsDeviceId() {
  /*! \brief Get the fastest device on the machine (with maximum GFLOPS).
   *
   * This function returns the identifier of the best available GPU (with
   * maximum GFLOPS)
   */
  int current_device = 0, cores_per_sm = 0;
  int max_compute_perf = 0, max_perf_device = 0;
  int device_count = 0, best_SM_arch = 0;
  cudaDeviceProp deviceProp;

  cudaGetDeviceCount(&device_count);

  // Find the best major SM Architecture GPU device
  while (current_device < device_count) {
    if (devices[current_device]->get_properties().major > best_SM_arch) {
      best_SM_arch = devices[current_device]->get_properties().major;
    }
    current_device++;
  }

  // Find the best CUDA capable GPU device
  current_device = device_count - 1;

  while (current_device >= 0) {
    deviceProp = devices[current_device]->get_properties();

    if (deviceProp.computeMode != cudaComputeModeProhibited) {
      if ((deviceProp.major == 9999) && (deviceProp.minor == 9999)) {
        cores_per_sm = 1;
      } else {
        cores_per_sm = ConvertSMVer2Cores(deviceProp.major, deviceProp.minor);
      }
      int compute_perf =
          deviceProp.multiProcessorCount * cores_per_sm * deviceProp.clockRate;

      if (compute_perf >= max_compute_perf) {
        // If we find GPU with SM major > 2, search only these
        if (best_SM_arch > 2) {
          // If our device==dest_SM_arch, choose this, or else pass
          if (deviceProp.major == best_SM_arch) {
            max_compute_perf = compute_perf;
            max_perf_device = current_device;
          }
        } else {
          max_compute_perf = compute_perf;
          max_perf_device = current_device;
        }
      }
    }
    --current_device;
  }
  return max_perf_device;
}

std::string carma_context::magma_info() {
  std::ostringstream stream;
#ifdef USE_MAGMA
  magma_int_t major, minor, micro;
  magma_version(&major, &minor, &micro);

  stream << "MAGMA " << (long long)major << "." << (long long)minor << "."
         << (long long)micro << ", " << (long long)(8 * sizeof(magma_int_t))
         << "-bit magma_int_t, " << (long long)(8 * sizeof(void *))
         << "-bit pointer.";
#else   // ifdef USE_MAGMA
  stream << "MAGMA not used";
#endif  // USE_MAGMA
  return stream.str();
}
