#include <carma_obj.h>

template<class T> CUDPPDatatype get_scandatatype()
{
  return CUDPP_INT;
}
template<> CUDPPDatatype get_scandatatype<float>()
{
  return CUDPP_FLOAT;
}

template<class T> CUDPPDatatype get_datatype()
{
  return CUDPP_UINT;
}
template<> CUDPPDatatype get_datatype<float>()
{
  return CUDPP_FLOAT;
}
template<> CUDPPDatatype get_datatype<int>()
{
  return CUDPP_INT;
}


template<class T> int carma_obj<T>::sort_init(bool keysOnly)
{
  this->keysOnly = keysOnly;

  //FIX ME ! this needs to be done on the GPU
  if ((!keysOnly) && (this->getValues() == NULL)){
    unsigned int *h_values = new unsigned int[this->nb_elem];
    for(int i = 0; i < this->nb_elem; i++) h_values[i] = i;
    if (this->getValues() == NULL) cudaMalloc((void **)&(this->values),this->nb_elem * sizeof(unsigned int));
    cudaMemcpy(this->values,h_values,this->nb_elem*sizeof(unsigned int),cudaMemcpyHostToDevice);
    delete h_values;
  }

  // initialize parameters based on present CUDA device
  CUDPPConfiguration config;
  config.algorithm = CUDPP_SORT_RADIX;
  config.datatype = get_datatype<T>();
  config.options = CUDPP_OPTION_KEY_VALUE_PAIRS;
  if (keysOnly) config.options = CUDPP_OPTION_KEYS_ONLY;

  CUDPPResult result = CUDPP_SUCCESS;  

  result = cudppPlan(&(this->mScanPlan), config, this->nb_elem, 1, 0);	

  if(result != CUDPP_SUCCESS) {
      printf("Error in sort plan creation\n");
  }

  return EXIT_SUCCESS;
}

template int caObjS::sort_init(bool keysOnly);
template int caObjI::sort_init(bool keysOnly);

template<class T> int carma_obj<T>::host2deviceInd(unsigned int *ind){
  cutilSafeCall(cudaMemcpy(this->values,ind,this->nb_elem * sizeof(unsigned int),cudaMemcpyHostToDevice));
  return EXIT_SUCCESS;
}

template int caObjS::host2deviceInd(unsigned int *ind);

template<class T> int carma_obj<T>::device2hostInd(unsigned int *ind){
  cutilSafeCall(cudaMemcpy(ind,this->values,this->nb_elem * sizeof(unsigned int),cudaMemcpyDeviceToHost));
  return EXIT_SUCCESS;
}

template int caObjS::device2hostInd(unsigned int *ind);

template<class T> int carma_obj<T>::device2deviceInd(unsigned int *src){
  if (this->getValues() == NULL) {
    cudaMalloc((void **)&(this->values),this->nb_elem * sizeof(unsigned int));
  }
  cutilSafeCall(cudaMemcpy(this->values,src,this->nb_elem * sizeof(unsigned int),cudaMemcpyDeviceToDevice));
  return EXIT_SUCCESS;
}

template int caObjS::device2deviceInd(unsigned int *src);

template<class T> int carma_obj<T>::sort()
{
  cudppSort(this->mScanPlan, this->d_data, (void*)this->values, 32, this->nb_elem);
  return EXIT_SUCCESS;
}

template int caObjS::sort();
template int caObjI::sort();

template<class T> int carma_obj<T>::config_scan(char *type, int dir, int incl)
{

  CUDPPConfiguration config;
  config.algorithm = CUDPP_SCAN;
  config.datatype = get_scandatatype<T>();
  config.op = CUDPP_ADD; // sum
  if (!strcmp(type, "max")) config.op = CUDPP_MAX;
  if (!strcmp(type, "min")) config.op = CUDPP_MIN;
  if (!strcmp(type, "mult")) config.op = CUDPP_MULTIPLY;
  if (!strcmp(type, "add")) config.op = CUDPP_ADD;

  CUDPPOption direction = CUDPP_OPTION_BACKWARD; // dir = -1
  CUDPPOption inclusivity = CUDPP_OPTION_INCLUSIVE; // incl = 1
  if (dir == -1) direction = CUDPP_OPTION_FORWARD;
  if (incl == 0) inclusivity = CUDPP_OPTION_EXCLUSIVE;
  config.options = direction | inclusivity;

  CUDPPResult result = CUDPP_SUCCESS;  

  result = cudppPlan(&(this->mScanPlan), config, this->nb_elem, 1, 0);	

  if(result != CUDPP_SUCCESS) {
    printf("Error in scan plan creation\n");
  }

  if (this->getOData() == NULL) {
    cutilSafeCall(cudaMalloc((void**)&(this->getOData()), sizeof(T)*this->nb_elem));
  }

  return EXIT_SUCCESS;
}

template int caObjS::config_scan(char *type, int dir, int incl);
template int caObjI::config_scan(char *type, int dir, int incl);

template<class T> T carma_obj<T>::scan()
{
  cudppScan(this->mScanPlan, this->o_data, this->d_data, this->nb_elem);

  T res;
  cutilSafeCallNoSync(cudaMemcpy((void **)&res,this->o_data,sizeof(T),cudaMemcpyDeviceToHost));
  return res;
}

template float caObjS::scan();
template int caObjI::scan();


template<class T> int carma_obj<T>::compact_init()
{
  if (this->getValues() == NULL){
    cudaMalloc((void **)&(this->values),this->nb_elem * sizeof(unsigned int));
  }

  // initialize parameters based on present CUDA device
  CUDPPConfiguration config;
  config.algorithm = CUDPP_COMPACT;
  config.datatype = get_datatype<T>();
  config.options = CUDPP_OPTION_FORWARD;

  CUDPPResult result = CUDPP_SUCCESS;  

  result = cudppPlan(&(this->mScanPlan), config, this->nb_elem, 1, 0);	

  if(result != CUDPP_SUCCESS)
    {
      printf("Error in plan creation\n");
    }

  return EXIT_SUCCESS;
}

template int caObjS::compact_init();
template int caObjI::compact_init();

template<class T> int carma_obj<T>::compact(carma_obj<T> *dest)
{
  cudppCompact(this->mScanPlan, dest->d_data, dest->d_numValid, this->d_data, this->values, this->nb_elem);

  return EXIT_SUCCESS;
}

template int caObjS::compact(caObjS *dest);
template int caObjI::compact(caObjI *dest);

template<class T> int carma_obj<T>::compact(carma_obj<T> *dest,unsigned int *values)
{
  cudppCompact(this->mScanPlan, dest->d_data, dest->d_numValid, this->d_data, values, this->nb_elem);

  return EXIT_SUCCESS;
}

template int caObjS::compact(caObjS *dest,unsigned int *values);
template int caObjI::compact(caObjI *dest,unsigned int *values);
