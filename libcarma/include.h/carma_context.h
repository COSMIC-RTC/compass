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
#include <memory>
#include <vector_types.h>
#include <cuda_runtime_api.h>

#include <carma_utils.h>
#include <carma_cublas.h>
#include <carma_cusparse.h>

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

	std::string getName() {
		return properties.name;
	}

	size_t getTotalMem() {
		return totalMem;
	}
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
private:
	int ndevice;
	std::vector<carma_device *> devices;
	int activeDevice;
	int** can_access_peer;

	/// singleton context
	static std::shared_ptr<carma_context> s_instance;

	carma_context();
	carma_context(int num_device);
	carma_context(int nb_devices, int32_t *devices);

	carma_context& operator=(const carma_context&) {
		return *s_instance;
	}
	carma_context(const carma_context&) :
			ndevice(0), activeDevice(0), can_access_peer(nullptr) {
	}

	void init_context(const int nb_devices, int32_t *devices_id);

public:
	~carma_context();

	static carma_context& instance_1gpu(int num_device);
	static carma_context& instance_ngpu(int nb_devices, int32_t *devices_id);
	static carma_context& instance();

	int get_ndevice() {
		return ndevice;
	}
	carma_device* get_device(int dev) {
		return devices[dev];
	}
	int get_activeDevice() {
		return activeDevice;
	}
	int get_activeRealDevice() {
		return devices[activeDevice]->get_id();
	}
	std::string get_DeviceName(int device);
	std::string get_DeviceInfo(int device);
	std::string get_DeviceMemInfo(int device);

	inline int _set_activeDeviceForCpy(int newDevice, int silent,
			std::string file, int line) {
		if (newDevice > ndevice)
			return -1;
		return (can_access_peer[activeDevice][newDevice] != 1) ?
				_set_activeDevice(newDevice, silent, file, line) : activeDevice;
	}
	inline int _set_activeDevice(int newDevice, int silent, std::string file,
			int line) {
		return (this->activeDevice != newDevice) ?
				_set_activeDeviceForce(newDevice, silent, file, line) :
				activeDevice;
	}
	int _set_activeDeviceForce(int newDevice, int silent, std::string file,
			int line);
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
	bool canP2P(int dev1, int dev2) {
		return can_access_peer[dev1][dev2];
	}

};

/// from /usr/local/cuda/samples/common/inc/helper_cuda.h
inline int ConvertSMVer2Cores(int major, int minor) {
	// Defines for GPU Architecture types (using the SM version to determine the # of cores per SM
	typedef struct {
		int SM; // 0xMm (hexidecimal notation), M = SM Major version, and m = SM minor version
		int Cores;
	} sSMtoCores;

	sSMtoCores nGpuArchCoresPerSM[] = { { 0x20, 32 }, // Fermi Generation (SM 2.0) GF100 class
			{ 0x21, 48 }, // Fermi Generation (SM 2.1) GF10x class
			{ 0x30, 192 }, // Kepler Generation (SM 3.0) GK10x class
			{ 0x32, 192 }, // Kepler Generation (SM 3.2) GK10x class
			{ 0x35, 192 }, // Kepler Generation (SM 3.5) GK11x class
			{ 0x37, 192 }, // Kepler Generation (SM 3.7) GK21x class
			{ 0x50, 128 }, // Maxwell Generation (SM 5.0) GM10x class
			{ 0x52, 128 }, // Maxwell Generation (SM 5.2) GM20x class
            { 0x53, 128}, // Maxwell Generation (SM 5.3) GM20x class
            { 0x60, 64 }, // Pascal Generation (SM 6.0) GP100 class
            { 0x61, 128}, // Pascal Generation (SM 6.1) GP10x class
            { 0x62, 128}, // Pascal Generation (SM 6.2) GP10x class
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
			major, minor, nGpuArchCoresPerSM[index - 1].Cores);
	return nGpuArchCoresPerSM[index - 1].Cores;
}
#endif // _CARMA_CONTEXT_H_
