#include <carma_context.h>
#include <iostream>
#include <sstream>
#include <sys/time.h>

#ifdef USE_CULA
// CULA headers
#include <cula.hpp>
#endif

#ifdef USE_MAGMA
// MAGMA headers
#include "magma.h"
#include "magma_lapack.h"
#endif

carma_device::carma_device(int devid) {
  this->id = devid;
  cudaGetDeviceProperties(&(this->properties), devid);
  this->sm_per_multiproc = ConvertSMVer2Cores(this->properties.major,
      this->properties.minor);
  this->compute_perf = this->properties.multiProcessorCount
      * this->sm_per_multiproc * this->properties.clockRate;

  this->p2p_activate = false;
}

carma_device::~carma_device() {
  this->id = -1;
}

carma_context::carma_context() {
  cutilSafeCall(cudaGetDeviceCount(&(this->ndevice)));
  if (this->ndevice == 0) {
    fprintf(stderr,
        "carma_context() CUDA error: no devices supporting CUDA.\n");
    throw "carma_context() CUDA error: no devices supporting CUDA.\n";
  }

  this->activeDevice = -1;
  int const size = ndevice;
  can_access_peer = new int*[size];
  for (int i = 0; i < ndevice; i++) {
    can_access_peer[i] = new int[size];
    for (int j = 0; j < ndevice; j++) {
      can_access_peer[i][j] = 0;
    }
  }

  int current_device = 0;
  int gpuid[64]; // we want to find the first two GPU's that can support P2P
  int gpu_count = 0; // GPUs that meet the criteria
  carma_device *current_yd = NULL;
  while (current_device < this->ndevice) {
    current_yd = new carma_device(current_device);
    devices.push_back(current_yd);

    if (current_yd->isGPUCapableP2P())
      gpuid[gpu_count++] = current_device;
    current_device++;
  }

  if (gpu_count > 1) {
    bool has_uva = true;
    for (int i = 0; i < gpu_count - 1; i++) {
      has_uva &= devices[gpuid[i]]->get_properties().unifiedAddressing;
      for (int j = i + 1; j < gpu_count; j++) {
        cutilSafeCall(
            cudaDeviceCanAccessPeer(&can_access_peer[gpuid[i]][gpuid[j]],
                gpuid[i], gpuid[j]));
        cutilSafeCall(
            cudaDeviceCanAccessPeer(&can_access_peer[gpuid[j]][gpuid[i]],
                gpuid[i], gpuid[j]));
        if ((can_access_peer[gpuid[i]][gpuid[j]] == 1)
            && (can_access_peer[gpuid[j]][gpuid[i]] == 1)) {
          printf("*** Enabling peer access between GPU%d and GPU%d... ***\n",
              gpuid[i], gpuid[j]);
          cutilSafeCall(cudaSetDevice(gpuid[i]));
          cutilSafeCall(cudaDeviceEnablePeerAccess(gpuid[j], 0));
          cutilSafeCall(cudaSetDevice(gpuid[j]));
          cutilSafeCall(cudaDeviceEnablePeerAccess(gpuid[i], 0));
        }
      }
    }
    has_uva &=
        devices[gpuid[gpu_count - 1]]->get_properties().unifiedAddressing;
    if (has_uva) {
      printf("*** All GPUs listed can support UVA... ***\n");
    }
  }

  this->activeDevice = set_activeDevice(get_maxGflopsDeviceId(), 1);

  carma_initCublas(&cublasHandle);
  carma_initCusparse(&cusparseHandle);

//  cusparsePointerMode_t mode;
//  cusparseGetPointerMode(cusparseHandle, &mode);
//  DEBUG_TRACE("%d\n", mode);

#ifdef USE_CULA
  // CULA init 
  culaStatus status = culaInitialize();
  if (status) {
    char buf[256];
    culaGetErrorInfoString(status, culaGetErrorInfo(), buf, sizeof(buf));
    printf("%s\n", buf);
  }
#endif

#ifdef USE_MAGMA
  // MAGMA init 
  magma_init();
#endif

#if DEBUG
  printf("CARMA Context created @ %8.8lX\n", (unsigned long)this);
#endif
}

carma_context::~carma_context() {
#ifdef USE_CULA
  // CULA finalize
  culaShutdown();
#endif

#ifdef USE_MAGMA
// MAGMA finalize
  magma_finalize();
#endif

  carma_shutdownCublas(cublasHandle);
  carma_shutdownCusparse(cusparseHandle);

  size_t idx = 0;
  while(this->devices.size()>0){
    delete this->devices.back();
    this->devices.pop_back();
    delete[] can_access_peer[idx++];
  }
  delete[] can_access_peer;
#if DEBUG
  printf("CARMA Context deleted @ %8.8lX\n", (unsigned long)this);
#endif
}

int carma_context::set_activeDeviceForCpy(int newDevice, int silent) {
  if (activeDevice == newDevice)
    return activeDevice;
  if (can_access_peer[activeDevice][newDevice] == 1)
    return activeDevice;
  return set_activeDevice(newDevice, silent);
}

int carma_context::set_activeDevice(int newDevice, int silent) {
  if (this->activeDevice == newDevice)
    return this->activeDevice;

  if (newDevice < ndevice) {
    cudaDeviceProp deviceProp;
    cutilSafeCall(cudaGetDeviceProperties(&deviceProp, newDevice));
    cutilSafeCall(cudaSetDevice(newDevice));
#ifdef USE_CULA
    culaStatus status = culaSelectDevice(newDevice);
    if(status) {
      char buf[256];
      culaGetErrorInfoString(status, culaGetErrorInfo(), buf, sizeof(buf));
      printf("%s\n", buf);
    }
#endif
    if (!silent)
      cout << "Using device " << newDevice << ": \"" << deviceProp.name
          << "\" with Compute " << deviceProp.major << "." << deviceProp.minor
          << " capability" << endl;
#if DEBUG
    cout << "DEBUG: Using device " << newDevice <<": \"" << deviceProp.name
    << "\" with Compute " << deviceProp.major << "."
    << deviceProp.minor << " capability" << endl;
#endif
    activeDevice = newDevice;
  } else {
    cerr << "Invalid Device Id : " << newDevice << " Your system has only "
        << ndevice << " CUDA capable device(s) available " << endl;
    cerr << "Leaving activeDevice to its current value : " << activeDevice
        << endl;
  }
  return activeDevice;
}

string carma_context::get_activeDeviceStr() {
  cudaDeviceProp deviceProp;
  cutilSafeCall(cudaGetDeviceProperties(&deviceProp, activeDevice));
  stringstream buf;
  buf << "Using device " << activeDevice << ": \"" << deviceProp.name
      << "\" with Compute " << deviceProp.major << "." << deviceProp.minor
      << " capability" << endl;
  return buf.str();
}

int carma_context::get_maxGflopsDeviceId()
/*! \brief Get the fastest device on the machine (with maximum GFLOPS).
 *
 * This function returns the identifier of the best available GPU (with maximum GFLOPS)
 */
{
  int current_device = 0, sm_per_multiproc = 0;
  int max_compute_perf = 0, max_perf_device = 0;
  int device_count = 0, best_SM_arch = 0;
  int arch_cores_sm[3] = { 1, 8, 32 };
  cudaDeviceProp deviceProp;

  cudaGetDeviceCount(&device_count);
  // Find the best major SM Architecture GPU device
  while (current_device < device_count) {
    cudaGetDeviceProperties(&deviceProp, current_device);
    if (deviceProp.major > best_SM_arch) {
      best_SM_arch = deviceProp.major;
    }
    current_device++;
  }

  // Find the best CUDA capable GPU device
  current_device = device_count - 1;
  while (current_device >= 0) {
    cudaGetDeviceProperties(&deviceProp, current_device);
    if (deviceProp.major == 9999 && deviceProp.minor == 9999) {
      sm_per_multiproc = 1;
    } else if (deviceProp.major <= 2) {
      sm_per_multiproc = arch_cores_sm[deviceProp.major];
    } else {
      sm_per_multiproc = arch_cores_sm[2];
    }

    int compute_perf = deviceProp.multiProcessorCount * sm_per_multiproc
        * deviceProp.clockRate;
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
    --current_device;
  }
  return max_perf_device;
}
