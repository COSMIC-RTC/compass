#include <carma_obj.h>
#include <cstdlib> /* required for randomize() and random() */

template <class T>
int carma_obj<T>::init_prng(long seed) {
  /*
    cudaDeviceProp deviceProperties =
    current_context->get_device(device)->get_properties(); const int maxThreads
    = deviceProperties.maxThreadsPerBlock; const int maxBlockDim =
    (deviceProperties.maxThreadsDim)[0]; int genPerBlock = std::min(maxThreads,
    maxBlockDim) / 2; int blockCount = deviceProperties.multiProcessorCount * 2;
  */
  int genPerBlock, blockCount;
  getNumBlocksAndThreads(current_context->get_device(device), nb_elem,
                         blockCount, genPerBlock);

  // Allocate memory for RNG states
  if (this->d_states == NULL)
    carmaSafeCall(cudaMalloc((void **)&(this->d_states),
                             blockCount * genPerBlock * sizeof(curandState)));
  carmaSafeCall(cudaMemset(this->d_states, 0,
                           blockCount * genPerBlock * sizeof(curandState)));

  this->nThreads = genPerBlock;
  this->nBlocks = blockCount;
  // randomize();
  std::vector<int> aseed(genPerBlock * blockCount);
  for (int cc = 0; cc <= genPerBlock * blockCount; cc++)
    aseed[cc] = seed + cc;  // random();

  int *seeds;
  carmaSafeCall(
      cudaMalloc((void **)&seeds, genPerBlock * blockCount * sizeof(int)));
  carmaSafeCall(cudaMemcpy(seeds, aseed.data(),
                           genPerBlock * blockCount * sizeof(int),
                           cudaMemcpyHostToDevice));

  // cerr << genPerBlock << " | " << blockCount << endl;
  carma_prng_init(seeds, genPerBlock, blockCount, this->d_states);

  carmaSafeCall(cudaFree(seeds));

  return EXIT_SUCCESS;
}

template int caObjI::init_prng(long seed);
template int caObjUI::init_prng(long seed);
template int caObjS::init_prng(long seed);
template int caObjD::init_prng(long seed);
template int caObjC::init_prng(long seed);
template int caObjZ::init_prng(long seed);

template <class T>
int carma_obj<T>::init_prng() {
  return this->init_prng(1234);
}

template int caObjI::init_prng();
template int caObjUI::init_prng();
template int caObjS::init_prng();
template int caObjD::init_prng();
template int caObjC::init_prng();
template int caObjZ::init_prng();

template <class T>
int carma_obj<T>::destroy_prng() {
  carmaSafeCall(cudaFree(this->d_states));
  return EXIT_SUCCESS;
}

template int caObjI::destroy_prng();
template int caObjUI::destroy_prng();
template int caObjS::destroy_prng();
template int caObjD::destroy_prng();
template int caObjC::destroy_prng();
template int caObjZ::destroy_prng();
template int carma_obj<tuple_t<float> >::destroy_prng();

template <class T>
int carma_obj<T>::prng(T *output, char gtype, float alpha, float beta) {
  carma_prng_cu(output, this->nThreads, this->nBlocks, this->d_states, gtype,
                this->nb_elem, alpha, beta);

  return EXIT_SUCCESS;
}

template int caObjI::prng(int *output, char gtype, float alpha, float beta);
template int caObjUI::prng(unsigned int *output, char gtype, float alpha,
                           float beta);
template int caObjS::prng(float *output, char gtype, float alpha, float beta);
template int caObjD::prng(double *output, char gtype, float alpha, float beta);
template int caObjC::prng(cuFloatComplex *output, char gtype, float alpha,
                          float beta);
template int caObjZ::prng(cuDoubleComplex *output, char gtype, float alpha,
                          float beta);

template <class T>
int carma_obj<T>::prng_montagn(float init_montagn) {
  carma_curand_montagn(this->d_states, this->d_data, this->nb_elem,
                       this->current_context->get_device(this->device));

  return EXIT_SUCCESS;
}

template int caObjI::prng_montagn(float init_montagn);
template int caObjUI::prng_montagn(float init_montagn);
template int caObjS::prng_montagn(float init_montagn);
template int caObjD::prng_montagn(float init_montagn);

template <class T>
int carma_obj<T>::prng(T *output, char gtype, float alpha) {
  carma_prng_cu(output, this->nThreads, this->nBlocks, this->d_states, gtype,
                this->nb_elem, alpha, 0.0f);

  return EXIT_SUCCESS;
}

template int caObjI::prng(int *output, char gtype, float alpha);
template int caObjUI::prng(unsigned int *output, char gtype, float alpha);
template int caObjS::prng(float *output, char gtype, float alpha);
template int caObjD::prng(double *output, char gtype, float alpha);
template int caObjC::prng(cuFloatComplex *output, char gtype, float alpha);
template int caObjZ::prng(cuDoubleComplex *output, char gtype, float alpha);

template <class T>
int carma_obj<T>::prng(char gtype) {
  return prng(this->d_data, gtype, 1.0f, 0.0f);
}

template int caObjI::prng(char gtype);
template int caObjUI::prng(char gtype);
template int caObjS::prng(char gtype);
template int caObjD::prng(char gtype);
template int caObjC::prng(char gtype);
template int caObjZ::prng(char gtype);

template <class T>
int carma_obj<T>::prng(char gtype, float alpha) {
  return prng(this->d_data, gtype, alpha, 0.0f);
}

template int caObjI::prng(char gtype, float alpha);
template int caObjUI::prng(char gtype, float alpha);
template int caObjS::prng(char gtype, float alpha);
template int caObjD::prng(char gtype, float alpha);
template int caObjC::prng(char gtype, float alpha);
template int caObjZ::prng(char gtype, float alpha);

template <class T>
int carma_obj<T>::prng(char gtype, float alpha, float beta) {
  return prng(this->d_data, gtype, alpha, beta);
}

template int caObjI::prng(char gtype, float alpha, float beta);
template int caObjUI::prng(char gtype, float alpha, float beta);
template int caObjS::prng(char gtype, float alpha, float beta);
template int caObjD::prng(char gtype, float alpha, float beta);
template int caObjC::prng(char gtype, float alpha, float beta);
template int caObjZ::prng(char gtype, float alpha, float beta);

template <class T>
int carma_obj<T>::init_prng_host(int seed) {
  if (this->gen != NULL) curandDestroyGenerator(this->gen);

  curandCreateGenerator(&(this->gen), CURAND_RNG_PSEUDO_XORWOW);
  // CURAND_RNG_PSEUDO_MTGP32
  // CURAND_RNG_PSEUDO_XORWOW
  curandSetPseudoRandomGeneratorSeed(this->gen, seed);

  return EXIT_SUCCESS;
}

template int caObjI::init_prng_host(int seed);
template int caObjUI::init_prng_host(int seed);
template int caObjS::init_prng_host(int seed);
template int caObjD::init_prng_host(int seed);
template int caObjC::init_prng_host(int seed);
template int caObjZ::init_prng_host(int seed);

template <class T>
int carma_obj<T>::prng_host(char gtype) {
  DEBUG_TRACE("Not implemented");
  return EXIT_FAILURE;
}
template int caObjI::prng_host(char gtype);
template int caObjUI::prng_host(char gtype);

template <>
int caObjS::prng_host(char gtype) {
  if (gtype == 'U')
    curandGenerateUniform(this->gen, this->d_data, this->nb_elem);
  if (gtype == 'N')
    curandGenerateNormal(this->gen, this->d_data, this->nb_elem, 0.0f, 1.0f);
  return EXIT_SUCCESS;
}

template <>
int caObjD::prng_host(char gtype) {
  if (gtype == 'U')
    curandGenerateUniformDouble(this->gen, this->d_data, this->nb_elem);
  if (gtype == 'N')
    curandGenerateNormalDouble(this->gen, this->d_data, this->nb_elem, 0.0,
                               1.0);
  return EXIT_SUCCESS;
}

template <class T>
int carma_obj<T>::prng_host(char gtype, T stddev) {
  DEBUG_TRACE("Not implemented");
  return EXIT_FAILURE;
}
template int caObjI::prng_host(char gtype, int stddev);
template int caObjUI::prng_host(char gtype, unsigned int stddev);

template <>
int caObjS::prng_host(char gtype, float stddev) {
  if (gtype == 'U')
    curandGenerateUniform(this->gen, this->d_data, this->nb_elem);
  if (gtype == 'N')
    curandGenerateNormal(this->gen, this->d_data, this->nb_elem, 0.0f, stddev);
  return EXIT_SUCCESS;
}
template <>
int caObjD::prng_host(char gtype, double stddev) {
  if (gtype == 'U')
    curandGenerateUniformDouble(this->gen, this->d_data, this->nb_elem);
  if (gtype == 'N')
    curandGenerateNormalDouble(this->gen, this->d_data, this->nb_elem, 0.0,
                               stddev);
  return EXIT_SUCCESS;
}

template <class T>
int carma_obj<T>::prng_host(char gtype, T stddev, T alpha) {
  DEBUG_TRACE("Not implemented");
  return EXIT_FAILURE;
}
template int caObjI::prng_host(char gtype, int stddev, int alpha);
template int caObjUI::prng_host(char gtype, unsigned int stddev,
                                unsigned int alpha);

template <>
int caObjS::prng_host(char gtype, float stddev, float alpha) {
  float *tmp_data;

  carmaSafeCall(cudaMalloc((void **)&tmp_data, this->nb_elem * sizeof(float)));

  if (gtype == 'U') curandGenerateUniform(this->gen, tmp_data, this->nb_elem);
  if (gtype == 'N')
    curandGenerateNormal(this->gen, tmp_data, this->nb_elem, 0.0f, stddev);

  carma_axpy(current_context->get_cublasHandle(), this->nb_elem, alpha,
             tmp_data, 1, this->d_data, 1);

  carmaSafeCall(cudaFree(tmp_data));

  return EXIT_SUCCESS;
}

template <>
int caObjD::prng_host(char gtype, double stddev, double alpha) {
  double *tmp_data;

  carmaSafeCall(cudaMalloc((void **)&tmp_data, this->nb_elem * sizeof(double)));

  if (gtype == 'U')
    curandGenerateUniformDouble(this->gen, tmp_data, this->nb_elem);
  if (gtype == 'N')
    curandGenerateNormalDouble(this->gen, tmp_data, this->nb_elem, 0.0f,
                               stddev);

  carma_axpy(current_context->get_cublasHandle(), this->nb_elem, alpha,
             tmp_data, 1, this->d_data, 1);

  carmaSafeCall(cudaFree(tmp_data));

  return EXIT_SUCCESS;
}

template <class T>
int carma_obj<T>::destroy_prng_host() {
  curandDestroyGenerator(this->gen);
  return EXIT_SUCCESS;
}

template int caObjS::destroy_prng_host();
template int caObjD::destroy_prng_host();
template int caObjC::destroy_prng_host();
template int caObjZ::destroy_prng_host();
template int caObjI::destroy_prng_host();
template int caObjUI::destroy_prng_host();
template int carma_obj<tuple_t<float> >::destroy_prng_host();
