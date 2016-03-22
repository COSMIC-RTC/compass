/**
 * \class carma_context
 *
 * \ingroup libcarma
 *
 * \brief this class provides the context in which carma_obj are created
 *
 * \author $Author: dg, as $
 *
 * \version $Revision: 1.0 $
 *
 * \date $Date: 2011/01/28$
 *
 */

#ifndef _CARMA_CONTEXT_H_
#define _CARMA_CONTEXT_H_

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <vector_types.h>
#include <cuda_runtime_api.h>

#include <carma_utils.h>
#include <carma_cublas.h>
#include <carma_cusparse.h>

using namespace std;

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

public:
  carma_device(int devid);
  carma_device(const carma_device& device);
  ~carma_device();

  int get_id() {
    return id;
  }
  cudaDeviceProp get_properties() {
    return properties;
  }
  float get_compute_perf() {
    return compute_perf;
  }
  float get_cores_per_sm() {
    return cores_per_sm;
  }
  bool isGPUCapableP2P() {
    return (bool) (properties.major >= 2);
  }

  bool isP2P_active() {
    return p2p_activate;
  }

  string getName() {return properties.name;}

  size_t getTotalMem() {return totalMem;}
  size_t getFreeMem() {
    carmaSafeCall(cudaMemGetInfo(&freeMem, &totalMem));
    return freeMem;
  }

  cublasHandle_t get_cublasHandle() {
    return cublasHandle;
  }
  cusparseHandle_t get_cusparseHandle() {
    return cusparseHandle;
  }
};

#define set_activeDevice(newDevice, silent)       _set_activeDevice(newDevice, silent, __FILE__, __LINE__)
#define set_activeDeviceForce(newDevice, silent)  _set_activeDeviceForce(newDevice, silent, __FILE__, __LINE__)
#define set_activeDeviceForCpy(newDevice, silent) _set_activeDeviceForCpy(newDevice, silent, __FILE__, __LINE__)

class carma_context {
protected:
  int ndevice;
  vector<carma_device *> devices;
  int activeDevice;
  int** can_access_peer;

  /// singleton context
  static carma_context *s_instance;

  carma_context();
  carma_context(int num_device);
  carma_context(const carma_context& cntxt);
public:
  ~carma_context();

  static carma_context *instance_1gpu(int num_device)
  {
      if (!s_instance)
        s_instance = new carma_context(num_device);
      return s_instance;
  }

  static carma_context *instance()
  {
      if (!s_instance)
        s_instance = new carma_context();
      return s_instance;
  }

  void kill(){
    if (s_instance)
      delete s_instance;
  }

  int get_ndevice() {
    return ndevice;
  }
  carma_device* get_device(int dev) {
    return devices[dev];
  }
  int get_activeDevice() {
    return activeDevice;
  }
  string get_DeviceName(int device);
  string get_DeviceInfo(int device);
  string get_DeviceMemInfo(int device);

  inline int _set_activeDeviceForCpy(int newDevice, int silent, string file, int line) {
    if(newDevice>ndevice) return -1;
    return (can_access_peer[activeDevice][newDevice] != 1)?_set_activeDevice(newDevice, silent, file, line):activeDevice;
  }
  inline int _set_activeDevice(int newDevice, int silent, string file, int line) {
    return (this->activeDevice != newDevice)?_set_activeDeviceForce(newDevice, silent, file, line):activeDevice;
  }
  int _set_activeDeviceForce(int newDevice, int silent, string file, int line);
  int get_maxGflopsDeviceId();
  cublasHandle_t get_cublasHandle() {
    return get_cublasHandle(activeDevice);
  }
  cusparseHandle_t get_cusparseHandle() {
    return get_cusparseHandle(activeDevice);
  }

  cublasHandle_t get_cublasHandle(int device) {
    return devices[device]->get_cublasHandle();
  }
  cusparseHandle_t get_cusparseHandle(int device) {
    return devices[device]->get_cusparseHandle();
  }

};

/// from /usr/local/cuda/samples/common/inc/helper_cuda.h
inline int ConvertSMVer2Cores(int major, int minor) {
  // Defines for GPU Architecture types (using the SM version to determine the # of cores per SM
  typedef struct {
    int SM; // 0xMm (hexidecimal notation), M = SM Major version, and m = SM minor version
    int Cores;
  } sSMtoCores;

  sSMtoCores nGpuArchCoresPerSM[] = { { 0x10, 8 }, // Tesla Generation (SM 1.0) G80 class
      { 0x11, 8 }, // Tesla Generation (SM 1.1) G8x class
      { 0x12, 8 }, // Tesla Generation (SM 1.2) G9x class
      { 0x13, 8 }, // Tesla Generation (SM 1.3) GT200 class
      { 0x20, 32 }, // Fermi Generation (SM 2.0) GF100 class
      { 0x21, 48 }, // Fermi Generation (SM 2.1) GF10x class
      { 0x30, 192 }, // Kepler Generation (SM 3.0) GK10x class
      { 0x35, 192 }, // Kepler Generation (SM 3.5) GK11x class
      { 0x37, 192}, // Kepler Generation (SM 3.7) GK21x class
      { 0x50, 128}, // Maxwell Generation (SM 5.0) GM10x class
      { 0x52, 128}, // Maxwell Generation (SM 5.2) GM20x class
      { -1, -1 } };

  int index = 0;

  while (nGpuArchCoresPerSM[index].SM != -1) {
    if (nGpuArchCoresPerSM[index].SM == ((major << 4) + minor)) {
      return nGpuArchCoresPerSM[index].Cores;
    }

    index++;
  }

  // If we don't find the values, we default use the previous one to run properly
  printf(
      "MapSMtoCores for SM %d.%d is undefined.  Default to use %d Cores/SM\n",
      major, minor, nGpuArchCoresPerSM[index-1].Cores);
  return nGpuArchCoresPerSM[index-1].Cores;
}
#endif // _CARMA_CONTEXT_H_
