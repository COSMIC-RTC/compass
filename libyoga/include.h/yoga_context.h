/**
 * \class yoga_context
 *
 * \ingroup libyoga
 *
 * \brief this class provides the context in which yoga_obj are created
 *
 * \author $Author: dg, as $
 *
 * \version $Revision: 1.0 $
 *
 * \date $Date: 2011/01/28$
 *
 */

#ifndef _YOGA_CONTEXT_H_
#define _YOGA_CONTEXT_H_

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <driver_types.h>
#include <vector_types.h>
#include <cuda_runtime_api.h>

#include <yoga_utils.h>

#ifdef MPI
#include <mpi.h>
#endif

using namespace std;

class yoga_device {
 protected:
  int            id;
  cudaDeviceProp properties;
  float          compute_perf;
  float          sm_per_multiproc;
  bool		 p2p_activate;

 public:
  yoga_device(int devid);
  yoga_device(int devid, float sm, float perf, bool p2p);
  yoga_device(const yoga_device& device);
  ~yoga_device();

  int            get_id(){ return id; }
  cudaDeviceProp get_properties(){ return properties; }
  float          get_compute_perf(){ return compute_perf; }
  float          get_sm_per_multiproc(){ return sm_per_multiproc; }
  bool           isGPUCapableP2P()
  {
    return (bool)(properties.major >= 2);
  }

  bool           isP2P_active()
  {
    return p2p_activate;
  }
  bool           set_p2pactive(bool active)
  {
    p2p_activate = active;
    return active;
  }

};

class yoga_host {
 protected:
  int                    id;
  int                    ndevices;
  vector<yoga_device *>  devices;
  int 			 activeDevice;
  int**			 can_access_peer;
  string                 hostname;
  bool		         ismaster;

 public:
  yoga_host(int hostid, string hostname);
  yoga_host(const yoga_host& hsts);
#ifdef MPI
  yoga_host(int nslave, MPI_Comm intercom);
#endif
  ~yoga_host();

  int          get_id() {return id;}
  int          get_ndevice() {return ndevices;}
  yoga_device *get_device(int dev) {return devices[dev];}
  int          get_activeDevice() {return activeDevice;}
  string       get_activeDeviceStr();
  int          set_activeDevice(int newDevice, int silent=1);
  int          set_activeDeviceForCpy(int newDevice, int silent=1);
  string       get_hostname(){return hostname;}
  bool         is_master(){return ismaster;}
  bool         is_p2pactive(int dev){return devices[dev]->isP2P_active();}

  int          get_maxGflopsDeviceId();

#ifdef MPI
  int          mpi_send_host(MPI_Comm intercom);
#endif
};

class yoga_context {
 protected:
  int                    ndevice;
  int                    nhosts;
  vector<yoga_host *>    hosts;

 public:
  yoga_context();
  yoga_context(bool with_mpi,string filename);
  yoga_context(const yoga_context& cntxt);
  ~yoga_context();

  int          get_ndevices() {return hosts[0]->get_ndevice();}
  int          get_ndevices(int host) {return hosts[host]->get_ndevice();}
  int          get_nhosts() {return nhosts;}
  yoga_host   *get_host(int host) {return hosts[host];}
  yoga_device *get_device(int dev) {return hosts[0]->get_device(dev);}
  yoga_device *get_device(int host, int dev) {return hosts[host]->get_device(dev);}
  int          get_activeDevice() {return hosts[0]->get_activeDevice();}
  string       get_activeDeviceStr(){return hosts[0]->get_activeDeviceStr();};
  string       get_hostname(int host){return hosts[host]->get_hostname();};
  int          get_activeDevice(int host) {return hosts[host]->get_activeDevice();}
  string       get_activeDeviceStr(int host){return hosts[host]->get_activeDeviceStr();}
  int          set_activeDevice(int newDevice, int silent=1){return hosts[0]->set_activeDevice(newDevice,silent);}
  //int          set_activeDevice(int host, int newDevice, int silent=1){return hosts[host]->set_activeDevice(newDevice,silent);}
  int          set_activeDeviceForCpy(int newDevice, int silent=1){return hosts[0]->set_activeDeviceForCpy(newDevice,silent);}
  //int          set_activeDeviceForCpy(int host, int newDevice, int silent=1){return hosts[host]->set_activeDeviceForCpy(newDevice,silent);}

  int          get_maxGflopsDeviceId(){return hosts[0]->get_maxGflopsDeviceId();}
  int          get_maxGflopsDeviceId(int host){return hosts[host]->get_maxGflopsDeviceId();}

  bool         is_p2pactive(int host, int dev){return hosts[host]->is_p2pactive(dev);}
};

inline int ConvertSMVer2Cores(int major, int minor)
{
  // Defines for GPU Architecture types (using the SM version to determine the # of cores per SM
  typedef struct {
    int SM; // 0xMm (hexidecimal notation), M = SM Major version, and m = SM minor version
    int Cores;
  } sSMtoCores;

  sSMtoCores nGpuArchCoresPerSM[] = 
    { { 0x10,  8 }, // Tesla Generation (SM 1.0) G80 class
      { 0x11,  8 }, // Tesla Generation (SM 1.1) G8x class
      { 0x12,  8 }, // Tesla Generation (SM 1.2) G9x class
      { 0x13,  8 }, // Tesla Generation (SM 1.3) GT200 class
      { 0x20, 32 }, // Fermi Generation (SM 2.0) GF100 class
      { 0x21, 48 }, // Fermi Generation (SM 2.1) GF10x class
      {   -1, -1 }
    };

  int index = 0;
  while (nGpuArchCoresPerSM[index].SM != -1) {
    if (nGpuArchCoresPerSM[index].SM == ((major << 4) + minor) ) {
      return nGpuArchCoresPerSM[index].Cores;
    }
    index++;
  }
  fprintf(stderr,"none of the GPU(s) can be used");
  return -1;
}

#endif // _YOGA_CONTEXT_H_
