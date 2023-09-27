// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      carma_context.cpp
//! \ingroup   libcarma
//! \class     CarmaContext
//! \brief     this class provides the context in which CarmaObj are created
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include <carma_context.h>
#include "carma_cublas.h"
#include "carma_cusparse.h"
#include "carma_cusolver.h"

#ifdef USE_MAGMA
// MAGMA headers
#include "magma.h"
#include "magma_lapack.h"
#endif  // USE_MAGMA

CarmaDevice::CarmaDevice(int devid) {
  carma_safe_call(cudaSetDevice(devid));

  // Instruct CUDA to yield its thread when waiting for results from the device.
  // This can increase latency when waiting for the device, but can increase the
  // performance of CPU threads performing work in parallel with the device. see
  // also: cudaDeviceScheduleAuto, cudaDeviceScheduleSpin,
  // cudaDeviceScheduleYield
  // carma_safe_call(cudaSetDeviceFlags(cudaDeviceScheduleYield));

  this->id = devid;
  cudaGetDeviceProperties(&(this->properties), devid);
  this->cores_per_sm =
      convert_sm_version2cores(this->properties.major, this->properties.minor);
  this->compute_perf = this->properties.multiProcessorCount *
                       this->cores_per_sm * this->properties.clockRate;

  carma_safe_call(cudaMemGetInfo(&total_mem, &total_mem));

  carma_init_cublas(&cublas_handle);
  carma_init_cusparse(&cusparse_handle);
  carma_init_cusolver(&cusolver_handle);

  // cudaStreamCreate(&main_stream);
  //  cusparsePointerMode_t mode;
  //  cusparseGetPointerMode(cusparse_handle, &mode);
  //  DEBUG_TRACE("%d\n", mode);

  // DEBUG_TRACE("done\n");
}
int CarmaDevice::set_cublas_math_mode(bool tensor) {
  if (tensor)
    carma_checkCublasStatus(
        cublasSetMathMode(cublas_handle, CUBLAS_TENSOR_OP_MATH));
  else
    carma_checkCublasStatus(
        cublasSetMathMode(cublas_handle, CUBLAS_DEFAULT_MATH));

  return EXIT_SUCCESS;
}

CarmaDevice::~CarmaDevice() {
  carma_shutdown_cublas(cublas_handle);
  carma_shutdown_cusparse(cusparse_handle);
  carma_shutdown_cusolver(cusolver_handle);
  // cudaStreamDestroy(main_stream);

  this->id = -1;
}

std::shared_ptr<CarmaContext> CarmaContext::s_instance;

CarmaContext &CarmaContext::instance_1gpu(int num_device) {
  if (!CarmaContext::s_instance) {
    CarmaContext::s_instance =
        std::shared_ptr<CarmaContext>(new CarmaContext(num_device));
  }
  return *CarmaContext::s_instance;
}

CarmaContext &CarmaContext::instance_ngpu(int nb_devices,
                                            int32_t *devices_id) {
  if (!CarmaContext::s_instance) {
    CarmaContext::s_instance = std::shared_ptr<CarmaContext>(
        new CarmaContext(nb_devices, devices_id));
  }
  return *CarmaContext::s_instance;
}

CarmaContext &CarmaContext::instance() {
  if (!CarmaContext::s_instance) {
    CarmaContext::s_instance =
        std::shared_ptr<CarmaContext>(new CarmaContext());
  }
  return *CarmaContext::s_instance;
}

CarmaContext::CarmaContext(int num_device) {
  can_access_peer = nullptr;
  this->active_device = -1;
  this->ndevice = -1;

  int devices[1];
  devices[0] = num_device;
  init_context(1, devices);
}

CarmaContext::CarmaContext() {
  carma_safe_call(cudaGetDeviceCount(&(this->ndevice)));
  can_access_peer = nullptr;
  this->active_device = -1;

  if (this->ndevice == 0) {
    DEBUG_TRACE("CarmaContext() CUDA error: no devices supporting CUDA.");
    throw std::runtime_error(
        "CarmaContext() CUDA error: no devices supporting CUDA.");
  }

  int const size = this->ndevice;
  int32_t devices_id[size];

  for (int i = 0; i < size; ++i) devices_id[i] = i;
  init_context(this->ndevice, devices_id);
}

CarmaContext::CarmaContext(int nb_devices, int32_t *devices_id) {
  can_access_peer = nullptr;
  this->active_device = -1;
  this->ndevice = -1;

  init_context(nb_devices, devices_id);
}

void CarmaContext::init_context(const int nb_devices, int32_t *devices_id) {
  // TODO(all) : why seed is initialized here ?
  srandom(1234);
  this->active_device = -1;

  int n_cuda_devices = 0;

  carma_safe_call(cudaGetDeviceCount(&n_cuda_devices));

  if (!n_cuda_devices) {
    DEBUG_TRACE("CarmaContext() CUDA error: no devices supporting CUDA.");
    throw std::runtime_error(
        "CarmaContext() CUDA error: no devices supporting CUDA.");
  }

  if (nb_devices > n_cuda_devices) {
    DEBUG_TRACE(
        "CarmaContext() CUDA error: not enough devices supporting CUDA. ask "
        "%d, available %d",
        nb_devices, n_cuda_devices);
    DEBUG_TRACE("CarmaContext() will be initialized on GPU 0 only");
    this->ndevice = 1;
  } else
    this->ndevice = nb_devices;
  int current_device = 0;

  while (current_device < this->ndevice) {
    devices.push_back(new CarmaDevice(devices_id[current_device]));
    current_device++;
  }

  can_access_peer = new int *[this->ndevice];

  for (int i = 0; i < ndevice; i++) {
    can_access_peer[i] = new int[this->ndevice];

    for (int j = 0; j < ndevice; j++) {
      can_access_peer[i][j] = (i == j);
    }
  }

#ifdef USE_P2P_ACCESS

  int gpuid[this->ndevice];  // we want to find the first two GPU's that can
                             // support P2P
  int gpu_count = 0;         // GPUs that meet the criteria
  current_device = 0;

  while (current_device < this->ndevice) {
    if (devices[current_device]->is_gpu_capable_p2p())
      gpuid[gpu_count++] = current_device;
    current_device++;
  }

  if (gpu_count > 1) {
    for (int i = 0; i < gpu_count - 1; i++) {
      for (int j = i + 1; j < gpu_count; j++) {
        carma_safe_call(cudaDeviceCanAccessPeer(
            &can_access_peer[gpuid[i]][gpuid[j]], devices_id[gpuid[i]],
            devices_id[gpuid[j]]));
        carma_safe_call(cudaDeviceCanAccessPeer(
            &can_access_peer[gpuid[j]][gpuid[i]], devices_id[gpuid[j]],
            devices_id[gpuid[i]]));

        if ((can_access_peer[gpuid[i]][gpuid[j]] == 1) &&
            (can_access_peer[gpuid[j]][gpuid[i]] == 1)) {
          printf("*** Enabling peer access between GPU%d and GPU%d... ***\n",
                 devices_id[gpuid[i]], devices_id[gpuid[j]]);
          carma_safe_call(cudaSetDevice(devices_id[gpuid[i]]));
          carma_safe_call(cudaDeviceEnablePeerAccess(devices_id[gpuid[j]], 0));
          carma_safe_call(cudaSetDevice(devices_id[gpuid[j]]));
          carma_safe_call(cudaDeviceEnablePeerAccess(devices_id[gpuid[i]], 0));
        }
      }
    }
  }
#endif  // USE_P2P_ACCESS

  this->active_device =
      set_active_device_force(0, 1);  // get_max_gflops_device_id(), 1);

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

CarmaContext::~CarmaContext() {
  carma_safe_call(cudaDeviceSynchronize());

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

int CarmaContext::_set_active_device_force(int new_device, int silent,
                                          std::string file, int line) {
  if (new_device < ndevice) {
    carma_safe_call(cudaSetDevice(devices[new_device]->get_id()));
#if DEBUG
    silent = 0;
#endif  // DEBUG

    if (!silent) {
      std::cout << "Using device " << devices[new_device]->get_id() << ": \""
                << devices[new_device]->get_properties().name
                << "\" with Compute "
                << devices[new_device]->get_properties().major << "."
                << devices[new_device]->get_properties().minor << " capability"
                << std::endl;
    }
    active_device = new_device;
  } else {
    fprintf(stderr,
            "[%s:%d] Invalid Device Id : %d, Your system has only %d CUDA "
            "capable device(s) available ",
            file.c_str(), line, new_device, ndevice);
    std::cerr << "Leaving active_device to its current value : " << active_device
              << std::endl;
  }
  return active_device;
}

std::string CarmaContext::get_device_name(int device) {
  return devices[device]->get_properties().name;
}

std::string CarmaContext::get_device_info(int device) {
  std::stringstream buf;

  buf << "device " << device << ": \"" << devices[device]->get_properties().name
      << "\" with Compute " << devices[device]->get_properties().major << "."
      << devices[device]->get_properties().minor << " capability";
  return buf.str();
}

std::string CarmaContext::get_device_mem_info(int device) {
  std::stringstream buf;
  size_t total_mem = devices[device]->get_total_mem() / 1024 / 1024;
  size_t usedMem = total_mem - devices[device]->get_free_mem() / 1024 / 1024;

  buf << "device " << device << ": \"" << devices[device]->get_properties().name
      << "\" memory used " << usedMem << "MB / " << total_mem << "MB ("
      << usedMem * 100. / total_mem << "%)";
  return buf.str();
}

int CarmaContext::get_max_gflops_device_id() {
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
        cores_per_sm = convert_sm_version2cores(deviceProp.major, deviceProp.minor);
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

std::string CarmaContext::magma_info() {
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
