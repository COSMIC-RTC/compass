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

carma_context *carma_context::s_instance=NULL;


#ifdef USE_MAGMA
extern int g_magma_devices_cnt;
extern struct magma_device_info* g_magma_devices;
extern int g_init;
extern magma_queue_t* g_null_queues;
struct magma_device_info
{
    size_t memory;
    magma_int_t cuda_arch;
};
#endif

carma_device::carma_device(int devid) {
  carmaSafeCall(cudaSetDevice(devid));
  this->id = devid;
  cudaGetDeviceProperties(&(this->properties), devid);
  this->cores_per_sm = ConvertSMVer2Cores(this->properties.major,
      this->properties.minor);
  this->compute_perf = this->properties.multiProcessorCount
      * this->cores_per_sm * this->properties.clockRate;

  this->p2p_activate = false;

  carmaSafeCall(cudaMemGetInfo(&freeMem, &totalMem));

  carma_initCublas(&cublasHandle);
  carma_initCusparse(&cusparseHandle);

  //  cusparsePointerMode_t mode;
  //  cusparseGetPointerMode(cusparseHandle, &mode);
  //  DEBUG_TRACE("%d\n", mode);

  //DEBUG_TRACE("done\n");
}

carma_device::~carma_device() {
  carma_shutdownCublas(cublasHandle);
  carma_shutdownCusparse(cusparseHandle);

  this->id = -1;
}

carma_context::carma_context(int num_device){
  //TODO : why seed is initialized here ?
  srandom(1234);
  carmaSafeCall(cudaGetDeviceCount(&(this->ndevice)));
  can_access_peer=nullptr;

  if (this->ndevice == 0) {
    fprintf(stderr,
        "carma_context() CUDA error: no devices supporting CUDA.\n");
    throw "carma_context() CUDA error: no devices supporting CUDA.\n";
  }

  if (this->ndevice < num_device) {
    fprintf(stderr,
        "carma_context() CUDA error: not enought devices supporting CUDA.\n");
    throw "carma_context() CUDA error: not enought devices supporting CUDA.\n";
  }

  can_access_peer = new int*[1];
  can_access_peer[0] = new int[1];
  can_access_peer[0][0] = 1;

  this->activeDevice = 0;
  carma_device *current_yd = new carma_device(num_device);
  devices.push_back(current_yd);
  this->ndevice = 1;
  this->activeDevice = set_activeDeviceForce(0, 1); //get_maxGflopsDeviceId(), 1);

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
  // MAGMA custom init
  g_magma_devices_cnt=1;
  magma_malloc_cpu( (void**) &g_magma_devices, g_magma_devices_cnt * sizeof(struct magma_device_info) );
  cudaDeviceProp prop;
  cudaGetDeviceProperties( &prop, num_device );
  g_magma_devices[0].memory    = prop.totalGlobalMem;
  g_magma_devices[0].cuda_arch = prop.major*100 + prop.minor*10;
  magma_device_t cdev;
  magma_getdevice( &cdev );
  magma_malloc_cpu( (void**) &g_null_queues, g_magma_devices_cnt * sizeof(magma_queue_t) );
  magma_queue_create_from_cuda( num_device, NULL, NULL, NULL, &g_null_queues[0] );
  magma_setdevice( cdev );
  magmablasSetKernelStream( g_null_queues[0] );
  g_init += 1;  // increment (init - finalize) count
  #if DEBUG
  //  magma_print_environment();
  #endif
#endif

#if DEBUG
  printf("CARMA Context created @ %p\n", this);
#endif

}

carma_context::carma_context() {

  //TODO : why seed is initialized here ?
  srandom(1234);
  carmaSafeCall(cudaGetDeviceCount(&(this->ndevice)));
  can_access_peer=nullptr;

  if (this->ndevice == 0) {
    fprintf(stderr,
        "carma_context() CUDA error: no devices supporting CUDA.\n");
    throw "carma_context() CUDA error: no devices supporting CUDA.\n";
  }

  int const size = ndevice;
  can_access_peer = new int*[size];
  for (int i = 0; i < ndevice; i++) {
    can_access_peer[i] = new int[size];
    for (int j = 0; j < ndevice; j++) {
      can_access_peer[i][j] = 0;
    }
  }

  this->activeDevice = -1;
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
#ifdef USE_UVA
  if (gpu_count > 1) {
    bool has_uva = true;
    for (int i = 0; i < gpu_count - 1; i++) {
      has_uva &= devices[gpuid[i]]->get_properties().unifiedAddressing;
      for (int j = i + 1; j < gpu_count; j++) {
        carmaSafeCall(
            cudaDeviceCanAccessPeer(&can_access_peer[gpuid[i]][gpuid[j]],
                gpuid[i], gpuid[j]));
        carmaSafeCall(
            cudaDeviceCanAccessPeer(&can_access_peer[gpuid[j]][gpuid[i]],
                gpuid[i], gpuid[j]));
        if ((can_access_peer[gpuid[i]][gpuid[j]] == 1)
            && (can_access_peer[gpuid[j]][gpuid[i]] == 1)) {
          printf("*** Enabling peer access between GPU%d and GPU%d... ***\n",
              gpuid[i], gpuid[j]);
          carmaSafeCall(cudaSetDevice(gpuid[i]));
          carmaSafeCall(cudaDeviceEnablePeerAccess(gpuid[j], 0));
          carmaSafeCall(cudaSetDevice(gpuid[j]));
          carmaSafeCall(cudaDeviceEnablePeerAccess(gpuid[i], 0));
        }
      }
    }
    has_uva &=
        devices[gpuid[gpu_count - 1]]->get_properties().unifiedAddressing;
    if (has_uva) {
      printf("*** All GPUs listed can support UVA... ***\n");
    }
  }
#endif

  this->activeDevice = set_activeDeviceForce(0, 1); //get_maxGflopsDeviceId(), 1);

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
  #if DEBUG
  //  magma_print_environment();
  #endif
#endif

#if DEBUG
  printf("CARMA Context created @ %p\n", this);
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

  size_t idx = 0;
  while (this->devices.size() > 0) {
    delete this->devices.back();
    this->devices.pop_back();
    if(can_access_peer!=nullptr) delete[] can_access_peer[idx++];
  }
  delete[] can_access_peer;

  s_instance= NULL;
#if DEBUG
  printf("CARMA Context deleted @ %p\n", this);
#endif
}

int carma_context::_set_activeDeviceForce(int newDevice, int silent,
    string file, int line) {
  if (newDevice < ndevice) {
    carmaSafeCall(cudaSetDevice(devices[newDevice]->get_id()));
#ifdef USE_CULA
    culaStatus status = culaSelectDevice(newDevice);
    if(status) {
      char buf[256];
      culaGetErrorInfoString(status, culaGetErrorInfo(), buf, sizeof(buf));
      printf("%s\n", buf);
    }
#endif
#if DEBUG
    silent=0;
#endif
    if (!silent) {
      cout << "Using device " << devices[newDevice]->get_id() << ": \""
          << devices[newDevice]->get_properties().name << "\" with Compute "
          << devices[newDevice]->get_properties().major << "."
          << devices[newDevice]->get_properties().minor << " capability"
          << endl;
    }
    activeDevice = newDevice;
  } else {
    fprintf(stderr,
        "[%s:%d] Invalid Device Id : %d, Your system has only %d CUDA capable device(s) available ",
        file.c_str(), line, newDevice, ndevice);
    cerr << "Leaving activeDevice to its current value : " << activeDevice
        << endl;
  }
  return activeDevice;
}

string carma_context::get_DeviceName(int device) {
  stringstream buf;
  return devices[device]->get_properties().name;
}

string carma_context::get_DeviceInfo(int device) {
  stringstream buf;
  buf << "device " << device << ": \"" << devices[device]->get_properties().name
      << "\" with Compute " << devices[device]->get_properties().major << "."
      << devices[device]->get_properties().minor << " capability";
  return buf.str();
}

string carma_context::get_DeviceMemInfo(int device) {
  stringstream buf;
  size_t totalMem = devices[device]->getTotalMem()/1024/1024;
  size_t usedMem = totalMem - devices[device]->getFreeMem()/1024/1024;
  buf << "device " << device << ": \"" << devices[device]->get_properties().name
      << "\" memory used " << usedMem << "MB / "
      << totalMem << "MB (" << usedMem*100./totalMem << "%)";
  return buf.str();
}

int carma_context::get_maxGflopsDeviceId()
/*! \brief Get the fastest device on the machine (with maximum GFLOPS).
 *
 * This function returns the identifier of the best available GPU (with maximum GFLOPS)
 */
{
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
      if (deviceProp.major == 9999 && deviceProp.minor == 9999) {
        cores_per_sm = 1;
      } else {
        cores_per_sm = ConvertSMVer2Cores(deviceProp.major,
            deviceProp.minor);
      }
      int compute_perf = deviceProp.multiProcessorCount * cores_per_sm
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
    }
    --current_device;
  }
  return max_perf_device;
}
