#include <yoga_context.h>
#include <iostream>
#include <sstream>
#include <sys/time.h>
#include <cula.hpp>
#include <fstream>

yoga_device::yoga_device(int devid)
{
  this->id     = devid;
  cudaGetDeviceProperties(&(this->properties), devid);
  this->sm_per_multiproc = ConvertSMVer2Cores(this->properties.major,this->properties.minor);
  this->compute_perf  = this->properties.multiProcessorCount * this->sm_per_multiproc *
    this->properties.clockRate;
}

yoga_device::yoga_device(int devid, float sm, float perf, bool p2p)
{
  this->id               = devid;
  this->sm_per_multiproc = sm;
  this->compute_perf     = perf;
  this->p2p_activate     = p2p;
}

yoga_device::~yoga_device()
{
  this->id     = -1;
}

yoga_host::yoga_host(int hostid, string hostname)
{
  this->id       = hostid;
  this->hostname = hostname;

  cutilSafeCall(cudaGetDeviceCount( &(this->ndevices) ));
  if (this->ndevices == 0) {
    fprintf(stderr, "yoga_context() CUDA error: no devices supporting CUDA.\n");
    throw "yoga_context() CUDA error: no devices supporting CUDA.\n";
  }
  
  this->activeDevice = -1;
  int const size=ndevices;
  can_access_peer = new int*[size];
  for(int i=0; i<ndevices; i++){
    can_access_peer[i] = new int[size];
    for(int j=0; j<ndevices; j++){
      can_access_peer[i][j] = 0;
    }
  }
  
  int current_device = 0;
  int gpuid[64]; // we want to find the first two GPU's that can support P2P
  int gpu_count = 0;   // GPUs that meet the criteria
  yoga_device *current_yd = NULL;
  while ( current_device < this->ndevices ) {
    current_yd = new yoga_device(current_device);
    devices.push_back(current_yd);
    if(current_yd->isGPUCapableP2P()) gpuid[gpu_count++] = current_device;
    current_device++;
  }
  
  if(gpu_count>1){
    bool has_uva = true;
    for(int i=0; i<gpu_count-1; i++){
      has_uva &= devices[gpuid[i]]->get_properties().unifiedAddressing;
      for(int j=i+1; j<gpu_count; j++){
	cutilSafeCall(cudaDeviceCanAccessPeer(&can_access_peer[gpuid[i]][gpuid[j]], gpuid[i], gpuid[j]));
	cutilSafeCall(cudaDeviceCanAccessPeer(&can_access_peer[gpuid[j]][gpuid[i]], gpuid[i], gpuid[j]));
	if ((can_access_peer[gpuid[i]][gpuid[j]] == 1) && (can_access_peer[gpuid[j]][gpuid[i]] == 1)){
	  printf("*** Enabling peer access between GPU%d and GPU%d... ***\n", gpuid[i], gpuid[j]);
	  cutilSafeCall(cudaSetDevice(gpuid[i]));
	  cutilSafeCall(cudaDeviceEnablePeerAccess(gpuid[j], 0));
	  cutilSafeCall(cudaSetDevice(gpuid[j]));
	  cutilSafeCall(cudaDeviceEnablePeerAccess(gpuid[i], 0));
	  devices[gpuid[i]]->set_p2pactive(true);
	  devices[gpuid[j]]->set_p2pactive(true);
	}
      }
    }
    has_uva &= devices[gpuid[gpu_count-1]]->get_properties().unifiedAddressing;
    if (has_uva) {
      printf("*** All GPUs listed can support UVA... ***\n");
    }
  }
  this->activeDevice = set_activeDevice(get_maxGflopsDeviceId(), 1);
}
/*
  cudaGetDeviceProperties(&(this->properties), devid);
  this->sm_per_multiproc = ConvertSMVer2Cores(this->properties.major,this->properties.minor);
  this->compute_perf  = this->properties.multiProcessorCount * this->sm_per_multiproc *
  this->properties.clockRate;
*/

#ifdef MPI
yoga_host::yoga_host(int nslave, MPI_Comm intercom)
{
  int nchar_host;

  MPI_Recv(&id, 1, MPI_INT, nslave, 0, intercom,MPI_STATUS_IGNORE);
  MPI_Recv(&ndevices, 1, MPI_INT, nslave, 0, intercom,MPI_STATUS_IGNORE);
  MPI_Recv(&activeDevice, 1, MPI_INT, nslave, 0, intercom,MPI_STATUS_IGNORE);
  MPI_Recv(&nchar_host, 1, MPI_INT, nslave, 0, intercom,MPI_STATUS_IGNORE);
  MPI_Recv((char*)hostname.c_str(), nchar_host, MPI_CHAR, nslave, 0, intercom,MPI_STATUS_IGNORE);
  MPI_Recv(&ismaster, 1, MPI_C_BOOL, nslave, 0, intercom,MPI_STATUS_IGNORE);
  MPI_Recv(&can_access_peer, ndevices*ndevices, MPI_INT, nslave, 0, intercom,MPI_STATUS_IGNORE);

  int did;
  float perf,sm;
  bool p2p;
  yoga_device *current_yd = NULL;

  for (int cc=0;cc<ndevices;cc++) {
    MPI_Recv(&did, 1, MPI_INT, nslave, 0, intercom,MPI_STATUS_IGNORE);
    MPI_Recv(&perf, 1, MPI_FLOAT, nslave, 0, intercom,MPI_STATUS_IGNORE);
    MPI_Recv(&sm, 1, MPI_FLOAT, nslave, 0, intercom,MPI_STATUS_IGNORE);
    MPI_Recv(&p2p, 1, MPI_C_BOOL, nslave, 0, intercom,MPI_STATUS_IGNORE);
    current_yd = new yoga_device(did,perf,sm,p2p);
    devices.push_back(current_yd);
  }
}
#endif

yoga_host::~yoga_host()
{
  this->id     = -1;
  for (size_t idx = 0; idx < (this->devices).size(); idx++) {
    devices.pop_back();
    delete[] can_access_peer[idx];
  }
  delete[] can_access_peer;
}

#ifdef MPI
int yoga_host::mpi_send_host(MPI_Comm intercom)
{
  int id = get_id();
  MPI_Send(&id, 1, MPI_INT, 0, 0, intercom);
  int ndev = get_ndevice();
  MPI_Send(&ndev, 1, MPI_INT, 0, 0, intercom);
  int actd = get_activeDevice();
  MPI_Send(&actd, 1, MPI_INT, 0, 0, intercom);
  int leng = hostname.length();
  MPI_Send(&leng, 1, MPI_INT, 0, 0, intercom);
  MPI_Send((char*)hostname.c_str(), hostname.length(), MPI_CHAR, 0, 0, intercom);
  bool mast = is_master();
  MPI_Send(&mast, 1 , MPI_C_BOOL, 0, 0, intercom);
  //int *canaccess = can_access_peer;
  MPI_Send(&can_access_peer, ndevices*ndevices , MPI_INT, 0, 0, intercom);

  for (int cc=0;cc<ndevices;cc++) {
    int did = devices[cc]->get_id();
    MPI_Send(&did, 1, MPI_INT, 0, 0, intercom);
    float cperf = devices[cc]->get_compute_perf();
    MPI_Send(&cperf, 1, MPI_FLOAT, 0, 0, intercom);
    float sm = devices[cc]->get_sm_per_multiproc();
    MPI_Send(&sm, 1, MPI_FLOAT, 0, 0, intercom);
    bool act = devices[cc]->isP2P_active();
    MPI_Send(&act, 1 , MPI_C_BOOL, 0, 0, intercom);
  }

  return EXIT_SUCCESS;
}
#endif

int yoga_host::set_activeDeviceForCpy(int newDevice, int silent)
{
  if(activeDevice==newDevice) return activeDevice;
  if(can_access_peer[activeDevice][newDevice]==1) return activeDevice;
  return set_activeDevice(newDevice, silent);
}

int yoga_host::set_activeDevice(int newDevice, int silent)
{
  if(this->activeDevice==newDevice) return this->activeDevice;

  if (newDevice < ndevices) {
    cudaDeviceProp deviceProp;
    cutilSafeCall( cudaGetDeviceProperties(&deviceProp, newDevice) );
    cutilSafeCall( cudaSetDevice(newDevice) );
    culaStatus status = culaSelectDevice(newDevice);
    if(status){
      char buf[256];
      culaGetErrorInfoString(status, culaGetErrorInfo(), buf, sizeof(buf));
      printf("%s\n", buf);
    }
    if(!silent)
      cout << "Using device " << newDevice <<": \"" << deviceProp.name
	   << "\" with Compute " << deviceProp.major << "."
	   << deviceProp.minor << " capability" << endl;
    activeDevice=newDevice;
  } else {
    cout << "Invalid Device Id : " << newDevice << " Your system has only "
	 << ndevices << " CUDA capable device(s) available " << endl ;
    cout << "Leaving activeDevice to its current value : " << activeDevice << endl;
  }
  return activeDevice;
}

string yoga_host::get_activeDeviceStr()
{
  cudaDeviceProp deviceProp;
  cutilSafeCall( cudaGetDeviceProperties(&deviceProp, activeDevice) );
  stringstream buf;
  buf << "Using device " << activeDevice <<": \"" << deviceProp.name
      << "\" with Compute " << deviceProp.major << "."
      << deviceProp.minor << " capability" << endl;
  return buf.str();
}

int yoga_host::get_maxGflopsDeviceId()
/*! \brief Get the fastest device on the machine (with maximum GFLOPS).
 *
 * This function returns the identifier of the best available GPU (with maximum GFLOPS)
 */
{
  int current_device   = 0, sm_per_multiproc = 0;
  int max_compute_perf = 0, max_perf_device  = 0;
  int device_count     = 0, best_SM_arch     = 0;
  int arch_cores_sm[3] = { 1, 8, 32 };
  cudaDeviceProp deviceProp;

  cudaGetDeviceCount( &device_count );
  // Find the best major SM Architecture GPU device
  while ( current_device < device_count ) {
    cudaGetDeviceProperties( &deviceProp, current_device );
    if (deviceProp.major > best_SM_arch) {
      best_SM_arch = deviceProp.major;
    }
    current_device++;
  }

  // Find the best CUDA capable GPU device
  current_device = device_count-1;
  while( current_device >=  0) {
    cudaGetDeviceProperties( &deviceProp, current_device );
    if (deviceProp.major == 9999 && deviceProp.minor == 9999) {
      sm_per_multiproc = 1;
    } else if (deviceProp.major <= 2) {
      sm_per_multiproc = arch_cores_sm[deviceProp.major];
    } else {
      sm_per_multiproc = arch_cores_sm[2];
    }

    int compute_perf  = deviceProp.multiProcessorCount * sm_per_multiproc *
      deviceProp.clockRate;
    if( compute_perf  >= max_compute_perf ) {
      // If we find GPU with SM major > 2, search only these
      if ( best_SM_arch > 2 ) {
	// If our device==dest_SM_arch, choose this, or else pass
	if (deviceProp.major == best_SM_arch) {
	  max_compute_perf  = compute_perf;
	  max_perf_device   = current_device;
	}
      } else {
	max_compute_perf  = compute_perf;
	max_perf_device   = current_device;
      }
    }
    --current_device;
  }
  return max_perf_device;
}

yoga_context::yoga_context()
{
  nhosts = 1;
  yoga_host *current_host = new yoga_host(0,"master");
  hosts.push_back(current_host);
}

yoga_context::yoga_context(bool with_mpi,string filename)
{
  if (with_mpi) {
#ifdef MPI
    int numprocs;
    int err=0, errcodes[256], rank;
    int namelen;

    MPI_Info info = MPI_INFO_NULL;
    
    MPI_Comm intercom;

    int nslaves=0;
    string hostname=" ";
    ifstream infile;

    MPI_Init(NULL,NULL);

    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    
    //MPI_Comm icom;

    infile.open ((char*)filename.c_str());
    while(hostname != "") {
      getline(infile,hostname);
      if (hostname != "") {
	nslaves += 1;  
      } 
    }
    infile.close();

    char foo2[] = "slave_context";
    //for (int i=0;i<nslaves;i++)
    err = MPI_Comm_spawn(foo2, MPI_ARGV_NULL, nslaves,
			 info, 0, MPI_COMM_SELF,
			 &intercom, errcodes);  
    if (err) printf("Error in MPI_Comm_spawn\n");
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Get_processor_name((char*)hostname.c_str(), &namelen);
    
    yoga_host *current_host = new yoga_host(0,hostname);
    hosts.push_back(current_host);

    for (int cc=0;cc<nslaves;cc++) { 
      current_host = new yoga_host(cc,intercom);
      hosts.push_back(current_host);
    }
    nhosts = nslaves + 1;
#endif
  } else {
    nhosts = 1;
    yoga_host *current_host = new yoga_host(0,"master");
    hosts.push_back(current_host);
  }
  ndevice = hosts[0]->get_ndevice();


#ifdef DEBUG
  printf("YOGA Context created @ %8.8lX\n", (unsigned long)this);
#endif
}

yoga_context::~yoga_context()
{
  for (size_t idx = 0; idx < (this->hosts).size(); idx++) {
    hosts.pop_back();
  }
#ifdef DEBUG
	printf("YOGA Context deleted @ %8.8lX\n", (unsigned long)this);
#endif
}

