#include <carma_obj.h>
#include <cstdlib> /* required for randomize() and random() */

template<class T>
int carma_obj<T>::init_prng(long seed) {
  cudaDeviceProp deviceProperties = current_context->get_device(device)->get_properties();
  int maxThreads = deviceProperties.maxThreadsPerBlock;
  int maxBlockDim = (deviceProperties.maxThreadsDim)[0];
  int genPerBlock = min(maxThreads, maxBlockDim) / 2;
  int blockCount = deviceProperties.multiProcessorCount * 2;

  // Allocate memory for RNG states
  carmaSafeCall(
      cudaMalloc((void ** )&(this->d_states),
          blockCount * genPerBlock * sizeof(curandState)));

  this->nThreads = genPerBlock;
  this->nBlocks = blockCount;
  //randomize();
  int aseed[genPerBlock * blockCount];
  for (int cc = 0; cc <= genPerBlock * blockCount; cc++)
    aseed[cc] = seed + cc;//random();

  int *seeds;
  carmaSafeCall(
      cudaMalloc((void ** )&seeds, genPerBlock * blockCount * sizeof(int)));
  cudaMemcpy(seeds, aseed, genPerBlock * blockCount * sizeof(int),
      cudaMemcpyHostToDevice);

  //cerr << genPerBlock << " | " << blockCount << endl;
  carma_prng_init((int *) seeds, genPerBlock, blockCount, this->d_states);
  cudaFree(seeds);

  return EXIT_SUCCESS;
}

template int
caObjS::init_prng(long seed);
template int
caObjD::init_prng(long seed);
template int
caObjC::init_prng(long seed);
template int
caObjZ::init_prng(long seed);

template <class T>
int carma_obj<T>::init_prng() {
	return this->init_prng(1234);
}

template int
caObjS::init_prng();
template int
caObjD::init_prng();
template int
caObjC::init_prng();
template int
caObjZ::init_prng();

template<class T>
int carma_obj<T>::destroy_prng() {
  carmaSafeCall(cudaFree(this->d_states));
  return EXIT_SUCCESS;
}

template int
caObjS::destroy_prng();
template int
caObjD::destroy_prng();
template int
caObjC::destroy_prng();
template int
caObjZ::destroy_prng();
template int
caObjI::destroy_prng();
template int
caObjUI::destroy_prng();
template int
carma_obj< tuple_t<float> >::destroy_prng();

template<class T>
int carma_obj<T>::prng(T *output, char gtype, float alpha, float beta) {
  carma_prng_cu(output, this->nThreads, this->nBlocks, this->d_states, gtype,
      this->nb_elem, alpha, beta);

  return EXIT_SUCCESS;
}

template int
caObjS::prng(float *output, char gtype, float alpha, float beta);
template int
caObjD::prng(double *output, char gtype, float alpha, float beta);
template int
caObjC::prng(cuFloatComplex *output, char gtype, float alpha, float beta);
template int
caObjZ::prng(cuDoubleComplex *output, char gtype, float alpha, float beta);

template<class T>
int carma_obj<T>::prng(T *output, char gtype, float alpha) {
  carma_prng_cu(output, this->nThreads, this->nBlocks, this->d_states, gtype,
      this->nb_elem, alpha, 0.0f);

  return EXIT_SUCCESS;
}

template int
caObjS::prng(float *output, char gtype, float alpha);
template int
caObjD::prng(double *output, char gtype, float alpha);
template int
caObjC::prng(cuFloatComplex *output, char gtype, float alpha);
template int
caObjZ::prng(cuDoubleComplex *output, char gtype, float alpha);

template<class T>
int carma_obj<T>::prng(char gtype) {
  return prng(this->d_data, gtype, 1.0f, 0.0f);
}

template int
caObjS::prng(char gtype);
template int
caObjD::prng(char gtype);
template int
caObjC::prng(char gtype);
template int
caObjZ::prng(char gtype);

template<class T>
int carma_obj<T>::prng(char gtype, float alpha) {
  return prng(this->d_data, gtype, alpha, 0.0f);
}

template int
caObjS::prng(char gtype, float alpha);
template int
caObjD::prng(char gtype, float alpha);
template int
caObjC::prng(char gtype, float alpha);
template int
caObjZ::prng(char gtype, float alpha);

template<class T>
int carma_obj<T>::prng(char gtype, float alpha, float beta) {
  return prng(this->d_data, gtype, alpha, beta);
}

template int
caObjS::prng(char gtype, float alpha, float beta);
template int
caObjD::prng(char gtype, float alpha, float beta);
template int
caObjC::prng(char gtype, float alpha, float beta);
template int
caObjZ::prng(char gtype, float alpha, float beta);

template<class T>
int carma_obj<T>::init_prng_host(int seed) {
  curandCreateGenerator(&(this->gen), CURAND_RNG_PSEUDO_MTGP32);
  //CURAND_RNG_PSEUDO_MTGP32
  //CURAND_RNG_PSEUDO_XORWOW
  curandSetPseudoRandomGeneratorSeed(this->gen, seed);

  return EXIT_SUCCESS;
}

template int
caObjS::init_prng_host(int seed);
template int
caObjD::init_prng_host(int seed);
template int
caObjC::init_prng_host(int seed);
template int
caObjZ::init_prng_host(int seed);

template<class T>
int carma_obj<T>::prng_host(char gtype) {
  return EXIT_FAILURE;
}

template<>
int caObjS::prng_host(char gtype) {
  if (gtype == 'U')
    curandGenerateUniform(this->gen, this->d_data, this->nb_elem);
  if (gtype == 'N')
    curandGenerateNormal(this->gen, this->d_data, this->nb_elem, 0.0f, 1.0f);
  return EXIT_SUCCESS;
}

template<>
int caObjD::prng_host(char gtype) {
  if (gtype == 'U')
    curandGenerateUniformDouble(this->gen, this->d_data, this->nb_elem);
  if (gtype == 'N')
    curandGenerateNormalDouble(this->gen, this->d_data, this->nb_elem, 0.0,
        1.0);
  return EXIT_SUCCESS;
}

template<class T>
int carma_obj<T>::prng_host(char gtype, T alpha) {
  return EXIT_FAILURE;
}

template<>
int caObjS::prng_host(char gtype, float alpha) {
  if (gtype == 'U')
    curandGenerateUniform(this->gen, this->d_data, this->nb_elem);
  if (gtype == 'N')
    curandGenerateNormal(this->gen, this->d_data, this->nb_elem, 0.0f, alpha);
  return EXIT_SUCCESS;
}
template<>
int caObjD::prng_host(char gtype, double alpha) {
  if (gtype == 'U')
    curandGenerateUniformDouble(this->gen, this->d_data, this->nb_elem);
  if (gtype == 'N')
    curandGenerateNormalDouble(this->gen, this->d_data, this->nb_elem, 0.0,
        alpha);
  return EXIT_SUCCESS;
}

template<class T>
int carma_obj<T>::destroy_prng_host() {
  curandDestroyGenerator(this->gen);
  return EXIT_SUCCESS;
}

template int
caObjS::destroy_prng_host();
template int
caObjD::destroy_prng_host();
template int
caObjC::destroy_prng_host();
template int
caObjZ::destroy_prng_host();
template int
caObjI::destroy_prng_host();
template int
caObjUI::destroy_prng_host();
template int
carma_obj< tuple_t <float> >::destroy_prng_host();

