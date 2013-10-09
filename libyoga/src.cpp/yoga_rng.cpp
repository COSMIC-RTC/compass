#include <yoga_obj.h>
#include <cstdlib> /* required for randomize() and random() */

template<class T> int yoga_obj<T>::init_prng(int device){
  struct cudaDeviceProp deviceProperties;
  // Get device properties
  cutilSafeCall(cudaGetDeviceProperties(&deviceProperties, device));
  
  int maxThreads  = deviceProperties.maxThreadsPerBlock;
  int maxBlockDim = (deviceProperties.maxThreadsDim)[0];
  int genPerBlock = min(maxThreads,maxBlockDim)/2;
  int blockCount  = deviceProperties.multiProcessorCount*2;

  // Allocate memory for RNG states
  cutilSafeCall(cudaMalloc((void **)&(this->d_states),blockCount * genPerBlock * 
			   sizeof(curandState)));

  this->nThreads = genPerBlock;
  this->nBlocks  = blockCount;
  //randomize();
  int aseed[genPerBlock*blockCount];
  for (int cc = 0;cc<=genPerBlock*blockCount;cc++) aseed[cc]= random();

  int *seeds;
  cutilSafeCall(cudaMalloc((void **)&seeds,genPerBlock*blockCount * sizeof(int)));
  cudaMemcpy(seeds,aseed,genPerBlock*blockCount * sizeof(int),cudaMemcpyHostToDevice); 

  //cerr << genPerBlock << " | " << blockCount << endl;
  yoga_prng_init((int *)seeds, genPerBlock, blockCount, this->d_states);
  cudaFree(seeds);

  return EXIT_SUCCESS;
}

template int yObjS::init_prng(int device);
template int yObjD::init_prng(int device);
template int yObjC::init_prng(int device);
template int yObjZ::init_prng(int device);

template<class T> int yoga_obj<T>::destroy_prng(){
  cutilSafeThreadSync();
  cutilSafeCall( cudaFree(this->d_states) ); 
  return EXIT_SUCCESS;
}

template int yObjS::destroy_prng();
template int yObjD::destroy_prng();
template int yObjC::destroy_prng();
template int yObjZ::destroy_prng();

template<class T> int yoga_obj<T>::prng(T *output,char gtype, float alpha, float beta)
{  
  yoga_prng_cu(output,this->nThreads,this->nBlocks,this->d_states,gtype,this->nb_elem,alpha,beta);

  return EXIT_SUCCESS;
}

template int yObjS::prng(float *output,char gtype,float alpha, float beta);
template int yObjD::prng(double *output,char gtype, float alpha, float beta);
template int yObjC::prng(cuFloatComplex *output,char gtype,float alpha, float beta);
template int yObjZ::prng(cuDoubleComplex *output,char gtype,float alpha, float beta);

template<class T> int yoga_obj<T>::prng(T *output,char gtype, float alpha)
{  
  yoga_prng_cu(output,this->nThreads,this->nBlocks,this->d_states,gtype,this->nb_elem,alpha,0.0f);

  return EXIT_SUCCESS;
}

template int yObjS::prng(float *output,char gtype,float alpha);
template int yObjD::prng(double *output,char gtype, float alpha);
template int yObjC::prng(cuFloatComplex *output,char gtype,float alpha);
template int yObjZ::prng(cuDoubleComplex *output,char gtype,float alpha);

template<class T> int yoga_obj<T>::prng(char gtype){
  return prng(this->d_data,gtype,1.0f,0.0f);
}

template int yObjS::prng(char gtype);
template int yObjD::prng(char gtype);
template int yObjC::prng(char gtype);
template int yObjZ::prng(char gtype);


template<class T> int yoga_obj<T>::prng(char gtype, float alpha){
  return prng(this->d_data,gtype,alpha,0.0f);
}

template int yObjS::prng(char gtype,float alpha);
template int yObjD::prng(char gtype, float alpha);
template int yObjC::prng(char gtype,float alpha);
template int yObjZ::prng(char gtype,float alpha);


template<class T> int yoga_obj<T>::prng(char gtype, float alpha,float beta){
  return prng(this->d_data,gtype,alpha,beta);
}

template int yObjS::prng(char gtype,float alpha,float beta);
template int yObjD::prng(char gtype,float alpha,float beta);
template int yObjC::prng(char gtype,float alpha,float beta);
template int yObjZ::prng(char gtype,float alpha,float beta);

template<class T> int yoga_obj<T>::init_prng_host(int seed)
{
  curandCreateGenerator(&(this->gen),CURAND_RNG_PSEUDO_MTGP32);
  //CURAND_RNG_PSEUDO_MTGP32
  //CURAND_RNG_PSEUDO_XORWOW
  curandSetPseudoRandomGeneratorSeed(this->gen ,seed);

  return EXIT_SUCCESS;
}

template int yObjS::init_prng_host(int seed);
template int yObjD::init_prng_host(int seed);
template int yObjC::init_prng_host(int seed);
template int yObjZ::init_prng_host(int seed);


template<class T> int yoga_obj<T>::prng_host(char gtype){
  return EXIT_FAILURE;
}

template<> int yObjS::prng_host(char gtype){
   if (gtype == 'U') 
     curandGenerateUniform(this->gen ,this-> d_data , this->nb_elem);
   if (gtype == 'N') 
     curandGenerateNormal(this->gen ,this-> d_data , this->nb_elem,0.0f,1.0f);
   return EXIT_SUCCESS;
}

template<> int yObjD::prng_host(char gtype){
   if (gtype == 'U') 
     curandGenerateUniformDouble(this->gen ,this-> d_data , this->nb_elem);
   if (gtype == 'N') 
     curandGenerateNormalDouble(this->gen ,this-> d_data , this->nb_elem,0.0,1.0);
   return EXIT_SUCCESS;
}

template<class T> int yoga_obj<T>::prng_host(char gtype, T alpha){
  return EXIT_FAILURE;
}

template<> int yObjS::prng_host(char gtype, float alpha){
   if (gtype == 'U') 
     curandGenerateUniform(this->gen ,this-> d_data , this->nb_elem);
   if (gtype == 'N')
     curandGenerateNormal(this->gen ,this-> d_data , this->nb_elem,0.0f,alpha);
   return EXIT_SUCCESS;
}
template<> int yObjD::prng_host(char gtype, double alpha){
   if (gtype == 'U') 
     curandGenerateUniformDouble(this->gen ,this-> d_data , this->nb_elem);
   if (gtype == 'N') 
     curandGenerateNormalDouble(this->gen ,this-> d_data , this->nb_elem,0.0,alpha);
   return EXIT_SUCCESS;
}

template<class T> int yoga_obj<T>::destroy_prng_host(){
  cutilSafeThreadSync();
  curandDestroyGenerator (this-> gen ) ; 
  return EXIT_SUCCESS;
}

template int yObjS::destroy_prng_host();
template int yObjD::destroy_prng_host();
template int yObjC::destroy_prng_host();
template int yObjZ::destroy_prng_host();
template int yObjI::destroy_prng_host();
template int yObjUI::destroy_prng_host();


