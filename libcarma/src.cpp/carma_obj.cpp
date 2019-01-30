#include <carma_obj.h>
#include <cstdlib> /* required for randomize() and random() */

template <class T_data>
carma_obj<T_data>::carma_obj(carma_context *current_context,
                             const long *dims_data) {
  /** \brief carma_obj creator.
   * \param current_context : the context in which the carma_obj is created
   * \param dims_data : the array size Yorick format : [ndims,dims1,dims2,...]
   */
  init(current_context, dims_data, NULL, true, 0);
}

template <class T_data>
carma_obj<T_data>::carma_obj(const carma_obj<T_data> *src) {
  /** \brief carma_obj creator.
   * \param src : carma_obj to copy
   */
  init(src->current_context, src->dims_data, src->d_data, false,
       src->get_nbStreams());
}

template <class T_data>
carma_obj<T_data>::carma_obj(carma_context *current_context,
                             const carma_obj<T_data> *src) {
  /** \brief carma_obj creator.
   * \param src : carma_obj to copy
   */
  init(current_context, src->dims_data, src->d_data, false, 0);
}

template <class T_data>
carma_obj<T_data>::carma_obj(carma_context *current_context,
                             const long *dims_data, const T_data *data) {
  /** \brief carma_obj creator.
   * \param dims_data : the array size Yorick format : [ndims,dims1,dims2,...]
   * \param data : the array
   */
  init(current_context, dims_data, data, true, 0);
}

template <class T_data>
carma_obj<T_data>::carma_obj(carma_context *current_context,
                             const long *dims_data, int nb_streams) {
  /** \brief carma_obj creator.
   * \param dims_data : the array size Yorick format : [ndims,dims1,dims2,...]
   */
  init(current_context, dims_data, NULL, true, nb_streams);
}

template <class T_data>
carma_obj<T_data>::carma_obj(carma_context *current_context,
                             const long *dims_data, const T_data *data,
                             int nb_streams) {
  /** \brief carma_obj creator.
   * \param dims_data : the array size Yorick format : [ndims,dims1,dims2,...]
   * \param data : the array
   */
  init(current_context, dims_data, data, true, nb_streams);
}

template <class T_data>
void carma_obj<T_data>::init(carma_context *context, const long *dims_data,
                             const T_data *data, bool fromHost,
                             int nb_streams) {
  this->current_context = context;
  const long size_data = dims_data[0] + 1;
  this->dims_data = new long[size_data];
  memcpy(this->dims_data, dims_data, size_data * sizeof(long));

  this->nb_elem = dims_data[1];
  for (int i = 2; i < size_data; i++) this->nb_elem *= dims_data[i];

  carmaSafeCall(
      cudaMalloc((void **)&(this->d_data), sizeof(T_data) * this->nb_elem));
  this->device = current_context->get_activeDevice();

  if (data == NULL)
    carmaSafeCall(cudaMemset(this->d_data, 0, sizeof(T_data) * this->nb_elem));

  else if (fromHost)
    this->host2device(data);
  else
    this->copyFrom(data, this->nb_elem);

  this->plan = 0;
  this->gen = 0;
  this->d_states = 0;
  this->values = 0;
  this->o_data = 0;
  carmaSafeCall(cudaMalloc((void **)&(this->d_numValid), sizeof(size_t)));

  this->streams = new carma_streams();
  this->add_stream(nb_streams);

#if DEBUG
  printf("CARMA Object created @ 0x%p on GPU%d\n", this,
         current_context->get_activeRealDevice());
#endif
}

/*
 template<class T_data>
 carma_obj<T_data>& carma_obj<T_data>::operator= (const carma_obj<T_data>& obj){

 if(this->d_data!=0) carmaSafeCall( cudaFree(this->d_data) );
 if(this->dims_data!=0) delete(this->dims_data);
 carma_obj<T_data> *new_obj = new carma_obj(carma_context *current_context,
 obj); return new_obj;
 }
 */

template <class T_data>
carma_obj<T_data>::~carma_obj() {
  /** \brief carma_obj destructor.
   */

  int old_device = current_context->get_activeDevice();
  current_context->set_activeDevice(this->device, 1);

  dealloc();
  this->d_data = 0;

  delete[](this->dims_data);
  this->dims_data = 0;

  cudaFree(this->d_numValid);
  this->d_numValid = 0;

  delete this->streams;

  if (this->values != 0) carmaSafeCall(cudaFree(this->values));

  // if (this->o_data!=0) carmaSafeCall( cudaFree(this->o_data) );

  /* Destroy the CUFFT plan. */
  if (this->plan != 0) carmafftSafeCall(cufftDestroy(this->plan));

  if (this->gen != 0) this->destroy_prng_host();

  if (this->d_states != 0) this->destroy_prng();
    //  carmaSafeCall(cudaFree(this->d_states));
    // this->d_states = 0;

#if DEBUG
  printf("CARMA Object deleted @ 0x%p on GPU%d\n", this,
         current_context->get_activeRealDevice());
#endif
  current_context->set_activeDevice(old_device, 1);
}

template <class T_data>
template <typename T_dest>
int carma_obj<T_data>::host2device(const T_dest *data) {
  /** \brief host2device data transfer.
   * \param data : input data
   *
   * this method fills d_input with the imput data
   */

  carmaSafeCall(cudaMemcpy(this->d_data, data, sizeof(T_dest) * this->nb_elem,
                           cudaMemcpyHostToDevice));

  return EXIT_SUCCESS;
}
template int carma_obj<unsigned int>::template host2device<unsigned int>(
    const unsigned int *data);
template int carma_obj<int>::template host2device<int>(const int *data);
template int carma_obj<float>::template host2device<float>(const float *data);
template int carma_obj<double>::template host2device<double>(
    const double *data);
template int carma_obj<cuFloatComplex>::template host2device<cuFloatComplex>(
    const cuFloatComplex *data);
template int carma_obj<cuDoubleComplex>::template host2device<cuDoubleComplex>(
    const cuDoubleComplex *data);

/*
 template<class T_data>
 T_data* carma_obj<T_data>::getData(){
 * \brief getData data transfer.
 * \return data : pointer on the data
 *

 return d_data;
 }

 template float* caObjS::getData();
 template double* caObjD::getData();
 template unsigned int* caObjUI::getData();
 template float2* caObjS2::getData();
 template double2* caObjD2::getData();
 */

template <class T_data>
template <typename T_dest>
int carma_obj<T_data>::device2host(T_dest *data) {
  /** \brief device2host data transfer.
   * \param data : output data
   *
   * this method copies the values in d_output to the output array
   */

  carmaSafeCall(cudaMemcpy(data, this->d_data, sizeof(T_dest) * this->nb_elem,
                           cudaMemcpyDeviceToHost));

  return EXIT_SUCCESS;
}

template int carma_obj<unsigned int>::template device2host<unsigned int>(
    unsigned int *data);
template int carma_obj<int>::template device2host<int>(int *data);
template int carma_obj<float>::template device2host<float>(float *data);
template int carma_obj<double>::template device2host<double>(double *data);
template int carma_obj<cuFloatComplex>::template device2host<cuFloatComplex>(
    cuFloatComplex *data);
template int carma_obj<cuDoubleComplex>::template device2host<cuDoubleComplex>(
    cuDoubleComplex *data);

#ifdef CAN_DO_HALF
template int carma_obj<half>::template host2device<half>(const half *data);
template int carma_obj<half>::template device2host<half>(half *data);
template <>
template <>
// std::enable_if_t<std::is_same<T_data, half>::value>
int carma_obj<half>::host2device<float>(const float *data) {
  copyFromFloatToHalf(data, this->d_data, this->nb_elem,
                      this->current_context->get_device(this->device));
}
template <>
template <>
// std::enable_if_t<std::is_same<T_data, half>::value>
int carma_obj<half>::device2host<float>(float *data) {
  copyFromHalfToFloat(this->d_data, data, this->nb_elem,
                      this->current_context->get_device(this->device));
}
#endif

template <class T_data>
int carma_obj<T_data>::host2deviceAsync(const T_data *data,
                                        cudaStream_t stream) {
  /** \brief host2device data transfer.
   * \param data : input data
   *
   * this method fills d_input with the imput data
   */

  carmaSafeCall(cudaMemcpyAsync(this->d_data, data,
                                sizeof(T_data) * this->nb_elem,
                                cudaMemcpyHostToDevice, stream));

  return EXIT_SUCCESS;
}

/*
 template<class T_data>
 T_data* carma_obj<T_data>::getData(){
 * \brief getData data transfer.
 * \return data : pointer on the data
 *

 return d_data;
 }

 template float* caObjS::getData();
 template double* caObjD::getData();
 template unsigned int* caObjUI::getData();
 template float2* caObjS2::getData();
 template double2* caObjD2::getData();
 */

template <class T_data>
int carma_obj<T_data>::device2hostAsync(T_data *data, cudaStream_t stream) {
  /** \brief device2host data transfer.
   * \param data : output data
   *
   * this method copies the values in d_output to the output array
   */
  cudaMemcpyAsync(data, this->d_data, sizeof(T_data) * this->nb_elem,
                  cudaMemcpyDeviceToHost, stream);

  return EXIT_SUCCESS;
}

template <class T_data>
int carma_obj<T_data>::device2hostOpt(T_data *data) {
  /** \brief device2host data transfer.
   * \param data : output data
   *
   * this method copies the values in d_output to the output array
   */
  if (this->o_data == 0) return EXIT_FAILURE;

  carmaSafeCall(cudaMemcpy(data, this->o_data, sizeof(T_data) * this->nb_elem,
                           cudaMemcpyDeviceToHost));

  return EXIT_SUCCESS;
}

template <class T_data>
int carma_obj<T_data>::copyInto(T_data *data, int nb_elem) {
  /** \brief device2host data transfer.
   * \param data : output data
   *
   * this method copies the values in d_output to the output array
   */
  if (nb_elem > this->nb_elem) nb_elem = this->nb_elem;

  carmaSafeCall(cudaMemcpy(data, this->d_data, sizeof(T_data) * nb_elem,
                           cudaMemcpyDeviceToDevice));

  return EXIT_SUCCESS;
}

template <class T_data>
int carma_obj<T_data>::copyFrom(const T_data *data, int nb_elem) {
  /** \brief device2host data transfer.
   * \param data : output data
   *
   * this method copies the values in d_output to the output array
   */
  if (nb_elem > this->nb_elem) nb_elem = this->nb_elem;

  carmaSafeCall(cudaMemcpy(this->d_data, data, sizeof(T_data) * nb_elem,
                           cudaMemcpyDeviceToDevice));

  return EXIT_SUCCESS;
}

template <class T_data>
T_data carma_obj<T_data>::sum() {
  return reduce<T_data>(this->d_data, this->nb_elem);
  // int nBlocks;
  // int nThreads;

  // this->current_context->set_activeDevice(device, 1);
  // sumGetNumBlocksAndThreads(this->nb_elem,
  //                           this->current_context->get_device(device),
  //                           nBlocks, nThreads);

  // reduce<T_data>(this->nb_elem, nThreads, nBlocks, this->d_data,
  // this->d_data);

  // carmaCheckMsg("Kernel execution failed");

  // sum partial block sums on GPU
  /*
  int s = nBlocks;
  while (s > 1) {
    int threads = 0, blocks = 0;
    this->current_context->set_activeDevice(device,1);
    sumGetNumBlocksAndThreads(s, this->current_context->get_device(device),
  blocks, threads);

    reduce<T_data>(s, threads, blocks, this->d_data, this->d_data);

    s = (s + (threads * 2 - 1)) / (threads * 2);

  }
  */
  // T_data *h_odata = new T_data[nBlocks];
  // cudaMemcpy(h_odata, this->d_data, sizeof(T_data), cudaMemcpyDeviceToHost);
  // return h_odata[0];
}

template <class T_data>
void carma_obj<T_data>::init_reduceCub() {
  init_reduceCubCU(this->cub_data, this->cub_data_size, this->d_data,
                   this->o_data, this->nb_elem);
}

template <class T_data>
void carma_obj<T_data>::reduceCub() {
  reduceCubCU(this->cub_data, this->cub_data_size, this->d_data, this->o_data,
              this->nb_elem);
}

template <class T_data>
void carma_obj<T_data>::clip(T_data min, T_data max) {
  clip_array<T_data>(this->d_data, min, max, this->nb_elem,
                     this->current_context->get_device(device));
}

template <class T_data>
int carma_obj<T_data>::transpose(carma_obj<T_data> *source) {
  transposeCU(source->d_data, this->d_data, this->dims_data[1],
              this->dims_data[2]);
  return EXIT_SUCCESS;
}

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

#ifdef CAN_DO_HALF
template <>
int caObjH::init_prng(long seed) {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_SUCCESS;
}
#endif

template <class T>
int carma_obj<T>::init_prng() {
  return this->init_prng(1234);
}

#ifdef CAN_DO_HALF
template <>
int caObjH::init_prng() {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_SUCCESS;
}
#endif

template <class T>
int carma_obj<T>::destroy_prng() {
  carmaSafeCall(cudaFree(this->d_states));
  return EXIT_SUCCESS;
}

#ifdef CAN_DO_HALF
template <>
int caObjH::destroy_prng() {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_SUCCESS;
}
#endif

template <class T>
int carma_obj<T>::prng(T *output, char gtype, float alpha, float beta) {
  carma_prng_cu(output, this->nThreads, this->nBlocks, this->d_states, gtype,
                this->nb_elem, alpha, beta);

  return EXIT_SUCCESS;
}

#ifdef CAN_DO_HALF
template <>
int caObjH::prng(half *output, char gtype, float alpha, float beta) {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_SUCCESS;
}
#endif

template <class T>
int carma_obj<T>::prng_montagn(float init_montagn) {
  carma_curand_montagn(this->d_states, this->d_data, this->nb_elem,
                       this->current_context->get_device(this->device));

  return EXIT_SUCCESS;
}

#ifdef CAN_DO_HALF
template <>
int caObjH::prng_montagn(float init_montagn) {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_SUCCESS;
}
#endif

template <class T>
int carma_obj<T>::prng(T *output, char gtype, float alpha) {
  carma_prng_cu(output, this->nThreads, this->nBlocks, this->d_states, gtype,
                this->nb_elem, alpha, 0.0f);

  return EXIT_SUCCESS;
}

#ifdef CAN_DO_HALF
template <>
int caObjH::prng(half *output, char gtype, float alpha) {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_SUCCESS;
}
#endif

template <class T>
int carma_obj<T>::prng(char gtype) {
  return prng(this->d_data, gtype, 1.0f, 0.0f);
}

#ifdef CAN_DO_HALF
template <>
int caObjH::prng(char gtype) {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_SUCCESS;
}
#endif

template <class T>
int carma_obj<T>::prng(char gtype, float alpha) {
  return prng(this->d_data, gtype, alpha, 0.0f);
}

#ifdef CAN_DO_HALF
template <>
int caObjH::prng(char gtype, float alpha) {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_SUCCESS;
}
#endif

template <class T>
int carma_obj<T>::prng(char gtype, float alpha, float beta) {
  return prng(this->d_data, gtype, alpha, beta);
}

#ifdef CAN_DO_HALF
template <>
int caObjH::prng(char gtype, float alpha, float beta) {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_SUCCESS;
}
#endif

template <class T>
int carma_obj<T>::init_prng_host(int seed) {
  if (this->gen != NULL) curandDestroyGenerator(this->gen);

  curandCreateGenerator(&(this->gen), CURAND_RNG_PSEUDO_XORWOW);
  // CURAND_RNG_PSEUDO_MTGP32
  // CURAND_RNG_PSEUDO_XORWOW
  curandSetPseudoRandomGeneratorSeed(this->gen, seed);

  return EXIT_SUCCESS;
}

#ifdef CAN_DO_HALF
template <>
int caObjH::init_prng_host(int seed) {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_SUCCESS;
}
#endif

template <class T>
int carma_obj<T>::prng_host(char gtype) {
  DEBUG_TRACE("Not implemented");
  return EXIT_FAILURE;
}

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

#ifdef CAN_DO_HALF
template <>
int caObjH::destroy_prng_host() {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_SUCCESS;
}
#endif

template <class T_data>
int carma_obj<T_data>::host2deviceVect(const T_data *data, int incx, int incy) {
  /** \brief host2device generic data transfer method for a blas object
   * vector. \param data : input data (x) \param incx : increment in x \param
   * incy : increment in y
   *
   * this method fills a vector with the input data
   */

  carma_checkCublasStatus(cublasSetVector(this->nb_elem, sizeof(T_data), data,
                                          incx, this->d_data, incy));
  return EXIT_SUCCESS;
}

template <class T_data>
int carma_obj<T_data>::device2hostVect(T_data *data, int incx, int incy) {
  /** \brief device2host generic data transfer method for a blas object
   * vector. \param data : output data (y) \param incx : increment in x \param
   * incy : increment in y
   *
   * this method fills output with the vector data
   */

  carma_checkCublasStatus(cublasGetVector(this->nb_elem, sizeof(T_data),
                                          this->d_data, incx, data, incy));
  return EXIT_SUCCESS;
}

template <class T_data>
int carma_obj<T_data>::host2deviceMat(const T_data *data, int lda, int ldb) {
  /** \brief host2device generic data transfer method.
   * \param data : input data  (A)
   * \param mat : matrix to fill(B)
   * \param lda : leading dim of A (# of rows)
   * \param ldb : leading dim of B (# of rows)
   *
   * this method fills mat with the input data
   */
  carma_checkCublasStatus(cublasSetMatrix(this->dims_data[1],
                                          this->dims_data[2], sizeof(T_data),
                                          data, lda, this->d_data, ldb));
  return EXIT_SUCCESS;
}

template <class T_data>
int carma_obj<T_data>::device2hostMat(T_data *data, int lda, int ldb) {
  /** \brief device2host generic data transfer.
   * \param data : output data  (A)
   * \param mat : matrix to copy(B)
   * \param lda : leading dim of A (# of rows)
   * \param ldb : leading dim of B (# of rows)
   *
   * this method fills output with the mat data
   */
  carma_checkCublasStatus(cublasGetMatrix(this->dims_data[1],
                                          this->dims_data[2], sizeof(T_data),
                                          this->d_data, lda, data, ldb));
  return EXIT_SUCCESS;
}

/*
 *  ____  _        _    ____  _
 * | __ )| |      / \  / ___|/ |
 * |  _ \| |     / _ \ \___ \| |
 * | |_) | |___ / ___ \ ___) | |
 * |____/|_____/_/   \_\____/|_|
 *
 */

template <class T_data>
int carma_obj<T_data>::aimax(int incx) {
  /** \brief aimax method
   * \param incx : increment in x
   *
   * this method finds the smallest index of the maximum magnitude element in
   * vect
   */
  return carma_where_amax(current_context->get_cublasHandle(), this->nb_elem,
                          this->d_data, incx);
}

template <class T_data>
int carma_obj<T_data>::aimin(int incx) {
  /** \brief aimin method
   * \param incx : increment in x
   *
   * this method finds the smallest index of the minimum magnitude element in
   * vect
   */
  return carma_where_amin(current_context->get_cublasHandle(), this->nb_elem,
                          this->d_data, incx);
};
template <class T_data>
T_data carma_obj<T_data>::asum(int incx) {
  /** \brief asum method
   * \param incx : increment in x
   *
   * this method computes the sum of the absolute values of the elements
   */
  return carma_getasum(current_context->get_cublasHandle(), this->nb_elem,
                       this->d_data, incx);
  ;
};
template <class T_data>
T_data carma_obj<T_data>::nrm2(int incx) {
  /** \brief getNrm2 method
   * \param n    : vect size
   * \param vect : vector to sum (x)
   * \param incx : increment in x
   *
   * this method computes the Euclidean norm of vect
   */
  return carma_nrm2(current_context->get_cublasHandle(), this->nb_elem,
                    this->d_data, incx);
}
template <class T_data>
T_data carma_obj<T_data>::dot(carma_obj<T_data> *source, int incx, int incy) {
  /** \brief ccomputeDot method
   * \param n    : vect size
   * \param vectx : first vector (x)
   * \param incx : increment in x
   * \param vecty : second vector (y)
   * \param incy : increment in y
   *
   * this method computes the dot product of two vectors
   */
  return carma_dot(current_context->get_cublasHandle(), this->nb_elem,
                   source->d_data, incx, this->d_data, incy);
}
template <class T_data>
void carma_obj<T_data>::scale(T_data alpha, int incx) {
  /** \brief vectScale method
   * \param n    : vect size
   * \param alpha : scale factor
   * \param vect : vector to scale (x)
   * \param incx : increment in x
   *
   * this method replaces vector x with alpha * x

   */
  carma_scal(current_context->get_cublasHandle(), this->nb_elem, alpha,
             this->d_data, incx);
}
template <class T_data>
void carma_obj<T_data>::swap(carma_obj<T_data> *source, int incx, int incy) {
  /** \brief vectSwap method
   * \param n    : vect size
   * \param vectx : first vector to swap (x)
   * \param incx : increment in x
   * \param vecty : second vector to swap (y)
   * \param incy : increment in y
   *
   * this method interchanges vector x with vector y
   */

  carma_swap(current_context->get_cublasHandle(), this->nb_elem, source->d_data,
             incx, this->d_data, incy);
}

template <class T_data>
void carma_obj<T_data>::copy(carma_obj<T_data> *source, int incx, int incy) {
  /** \brief vectCopy method
   * \param n    : vect size
   * \param vectx : vector to copy (x)
   * \param incx : increment in x
   * \param vecty : vector to fill (y)
   * \param incy : increment in y
   *
   * this method interchanges vector x with vector y
   */

  carma_copy(current_context->get_cublasHandle(), this->nb_elem, source->d_data,
             incx, this->d_data, incy);
};
template <class T_data>
void carma_obj<T_data>::axpy(T_data alpha, carma_obj<T_data> *source, int incx,
                             int incy) {
  /** \brief computeAXPY method
   * \param alpha : scale factor
   * \param vectx : first vector (x)
   * \param incx : increment in x
   * \param vecty : second vector (y)
   * \param incy : increment in y
   *
   * this method multiplies vector x by scalar alpha and adds the result to
   * vector y
   */
  carma_axpy(current_context->get_cublasHandle(), this->nb_elem, alpha,
             source->d_data, incx, this->d_data, incy);
}
template <class T_data>
void carma_obj<T_data>::rot(carma_obj<T_data> *source, int incx, int incy,
                            T_data sc, T_data ss) {
  /** \brief ccomputeRot method
   * \param n    : vect size
   * \param vectx : first vector (x)
   * \param incx : increment in x
   * \param vecty : second vector (y)
   * \param incx : increment in y
   * \param sc : cosinus of rotation angle
   * \param ss : sinus of rotation angle
   *
   * this method computes the dot product of two vectors
   */
  carma_rot(current_context->get_cublasHandle(), this->nb_elem, source->d_data,
            incx, this->d_data, incy, sc, ss);
}

/*
 *  ____  _        _    ____ ____
 * | __ )| |      / \  / ___|___ \
 * |  _ \| |     / _ \ \___ \ __) |
 * | |_) | |___ / ___ \ ___) / __/
 * |____/|_____/_/   \_\____/_____|
 *
 */

template <class T_data>
void carma_obj<T_data>::gemv(char trans, T_data alpha, carma_obj<T_data> *matA,
                             int lda, carma_obj<T_data> *vectx, int incx,
                             T_data beta, int incy) {
  /** \brief gemv method.
   * \param trans : type of op 'n' : nothing  / 't' or 'c' : transpose
   * \param alpha : alpha
   * \param matA : matrix A
   * \param lda : leading dim of A (# of rows)
   * \param vectx : vector x
   * \param incx : increment for x
   * \param incy : increment for y
   *
   * this method performs one of the matrix‐vector operations y = alpha *
   * op(A)
   * * x + beta * y
   */
  carma_gemv(current_context->get_cublasHandle(), trans, matA->dims_data[1],
             matA->dims_data[2], alpha, matA->d_data, lda, vectx->d_data, incx,
             beta, this->d_data, incy);
}
template <class T_data>
void carma_obj<T_data>::ger(T_data alpha, carma_obj<T_data> *vectx, int incx,
                            carma_obj<T_data> *vecty, int incy, int lda) {
  /** \brief ger method.
   * \param alpha : alpha
   * \param vectx : m-element vector
   * \param incx : increment for x
   * \param vecty : vector y
   * \param incy : increment for y
   * \param lda : leading dim of A (# of rows)
   *
   * this method performs the symmetric rank 1 operation A = alpha * x * y T +
   * A
   */
  carma_ger(current_context->get_cublasHandle(), this->dims_data[1],
            this->dims_data[2], alpha, vectx->d_data, incx, vecty->d_data, incy,
            this->d_data, lda);
}

template <class T_data>
void carma_obj<T_data>::symv(char uplo, T_data alpha, carma_obj<T_data> *matA,
                             int lda, carma_obj<T_data> *vectx, int incx,
                             T_data beta, int incy) {
  /** \brief symv method.
   * \param uplo : upper or lower fill mode
   * \param alpha : alpha
   * \param matA : matrix A
   * \param lda : leading dim of A (# of rows)
   * \param vectx : vector x
   * \param incx : increment for x
   * \param incy : increment for y
   *
   * this method performs one of the matrix‐vector operations y = alpha *
   * op(A)
   * * x + beta * y
   */
  carma_symv(current_context->get_cublasHandle(), uplo, matA->dims_data[1],
             alpha, matA->d_data, lda, vectx->d_data, incx, beta, this->d_data,
             incy);
}

/*
 *  ____  _        _    ____ _____
 * | __ )| |      / \  / ___|___ /
 * |  _ \| |     / _ \ \___ \ |_ \
 * | |_) | |___ / ___ \ ___) |__) |
 * |____/|_____/_/   \_\____/____/
 *
 */

template <class T_data>
void carma_obj<T_data>::gemm(char transa, char transb, T_data alpha,
                             carma_obj<T_data> *matA, int lda,
                             carma_obj<T_data> *matB, int ldb, T_data beta,
                             int ldc) {
  /** \brief generic gemm method.
   * \param transa : type of op 'n' : nothing  / 't' or 'c' : transpose
   * \param transb : type of op 'n' : nothing  / 't' or 'c' : transpose
   * \param alpha : alpha
   * \param matA : matrix A
   * \param lda : leading dim of A (# of rows)
   * \param matB : matrix B
   * \param ldb : leading dim of B (# of rows)
   * \param beta : beta
   * \param ldc : leading dim of C (# of rows)
   *
   * this method performs one of the matrix‐matrix operations:
   * C = alpha * op(A) * op(B) + beta * C ,
   * where  op(X) = X  or  op(X) = X^T
   */
  int k = (((transa == 'N') || (transa == 'n')) ? matA->dims_data[2]
                                                : matA->dims_data[1]);
  carma_gemm(current_context->get_cublasHandle(), transa, transb,
             this->dims_data[1], this->dims_data[2], k, alpha, matA->d_data,
             lda, matB->d_data, ldb, beta, this->d_data, ldc);
}
template <class T_data>
void carma_obj<T_data>::symm(char side, char uplo, T_data alpha,
                             carma_obj<T_data> *matA, int lda,
                             carma_obj<T_data> *matB, int ldb, T_data beta,
                             int ldc) {
  /** \brief generic symm method.
   * \param side : which side of the equation is symmetric matrix A
   * \param uplo : fill mode of matrix A
   * \param alpha : alpha
   * \param matA : matrix A
   * \param lda : leading dim of A (# of rows)
   * \param matB : matrix B
   * \param ldb : leading dim of B (# of rows)
   * \param beta : beta
   * \param ldc : leading dim of C (# of rows)
   *
   * this method performs one of the symmetric matrix‐matrix operations:
   * C = alpha * A * B + beta * C ,
   * or
   * C = alpha * B * A + beta * C ,
   * where A is a symmetric matrix
   */
  carma_symm(current_context->get_cublasHandle(), side, uplo,
             this->dims_data[1], this->dims_data[2], alpha, matA->d_data, lda,
             matB->d_data, ldb, beta, this->d_data, ldc);
}
template <class T_data>
void carma_obj<T_data>::syrk(char uplo, char transa, T_data alpha,
                             carma_obj<T_data> *matA, int lda, T_data beta,
                             int ldc) {
  /** \brief generic syrk method.
   * \param uplo : fill mode of matrix A
   * \param alpha : alpha
   * \param matA : matrix A
   * \param lda : leading dim of A (# of rows)
   * \param beta : beta
   * \param ldc : leading dim of C (# of rows)
   *
   * this method performs one of the symmetric matrix‐matrix operations:
   * C = alpha * op(A) * transpose(op(A)) + beta * C ,
   * where  op(X) = X  or  op(X) = X^T
   */
  /*
   carma_syrk(current_context->get_cublasHandle(), uplo, transa,
   this->dims_data[1], this->dims_data[2], alpha, matA->d_data, lda, beta,
   this->d_data, ldc);
   */
  carma_syrk(current_context->get_cublasHandle(), uplo, transa,
             matA->dims_data[1], matA->dims_data[2], alpha, matA->d_data, lda,
             beta, this->d_data, ldc);
  // "this" refer to matrix "C": this->dims_data[1]=this->dims_data[2]=n
}
template <class T_data>
void carma_obj<T_data>::syrkx(char uplo, char transa, T_data alpha,
                              carma_obj<T_data> *matA, int lda,
                              carma_obj<T_data> *matB, int ldb, T_data beta,
                              int ldc) {
  /** \brief generic syrkx method.
   * \param uplo : fill mode of matrix A
   * \param transa : type of op 'n' : nothing  / 't' or 'c' : transpose
   * \param alpha : alpha
   * \param matA : matrix A
   * \param lda : leading dim of A (# of rows)
   * \param matB : matrix B
   * \param ldb : leading dim of B (# of rows)
   * \param beta : beta
   * \param ldc : leading dim of C (# of rows)
   *
   * this method performs one of the symmetric matrix‐matrix operations:
   * C = alpha * op(A) * transpose(op(B)) + beta * C ,
   * where  op(X) = X  or  op(X) = X^T
   */
  /*carma_syrkx(current_context->get_cublasHandle(), uplo, transa,
   this->dims_data[1], this->dims_data[2], alpha, matA->d_data, lda,
   matB->d_data, ldb, beta, this->d_data, ldc);*/
  carma_syrkx(current_context->get_cublasHandle(), uplo, transa,
              matA->dims_data[1], matA->dims_data[2], alpha, matA->d_data, lda,
              matB->d_data, ldb, beta, this->d_data, ldc);
}
template <class T_data>
void carma_obj<T_data>::geam(char transa, char transb, T_data alpha,
                             carma_obj<T_data> *matA, int lda, T_data beta,
                             carma_obj<T_data> *matB, int ldb, int ldc) {
  /** \brief generic geam method.
   * \param transa : type of op 'n' : nothing  / 't' or 'c' : transpose
   * \param transb : type of op 'n' : nothing  / 't' or 'c' : transpose
   * \param alpha : alpha
   * \param matA : matrix A
   * \param lda : leading dim of A (# of rows)
   * \param beta : beta
   * \param matB : matrix B
   * \param ldb : leading dim of B (# of rows)
   * \param ldc : leading dim of C (# of rows)
   *
   * C = alpha * op(A) + beta * op(B),
   * where  op(X) = X  or  op(X) = X^T
   */

  carma_geam(current_context->get_cublasHandle(), transa, transb,
             this->dims_data[1], this->dims_data[2], alpha, matA->d_data, lda,
             beta, matB->d_data, ldb, this->d_data, ldc);
}
template <class T_data>
void carma_obj<T_data>::dgmm(char side, carma_obj<T_data> *matA, int lda,
                             carma_obj<T_data> *vectx, int incx, int ldc) {
  /** \brief generic dgmm method.
   * \param side : side of equation for matrix A
   * \param matA : matrix A
   * \param lda : leading dim of A (# of rows)
   * \param vectx : vector x
   * \param incx : increment on x
   * \param ldc : leading dim of C (# of rows)
   *
   * C = A * diag(X) or C = diag(X) * A
   */

  carma_dgmm(current_context->get_cublasHandle(), side, this->dims_data[1],
             this->dims_data[2], matA->d_data, lda, vectx->d_data, incx,
             this->d_data, ldc);
}

template class carma_obj<int>;
template class carma_obj<unsigned int>;
template class carma_obj<float>;
template class carma_obj<double>;
template class carma_obj<cuFloatComplex>;
template class carma_obj<cuDoubleComplex>;
// template class carma_obj<struct tuple_t<float>>;

#ifdef CAN_DO_HALF
template <>
int caObjH::aimax(int incx) {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_FAILURE;
}
template <>
int caObjH::aimin(int incx) {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_FAILURE;
}
template <>
half caObjH::asum(int incx) {
  DEBUG_TRACE("Not implemented for half precision");
  return __float2half(0.f);
}
template <>
half caObjH::nrm2(int incx) {
  DEBUG_TRACE("Not implemented for half precision");
  return __float2half(0.f);
}
template <>
void caObjH::scale(half alpha, int incx) {
  DEBUG_TRACE("Not implemented for half precision");
}
template <>
void caObjH::copy(caObjH *, int incx, int incy) {
  DEBUG_TRACE("Not implemented for half precision");
}
template <>
void caObjH::swap(caObjH *, int incx, int incy) {
  DEBUG_TRACE("Not implemented for half precision");
}
template <>
void caObjH::axpy(half alpha, caObjH *source, int incx, int incy) {
  DEBUG_TRACE("Not implemented for half precision");
}
template <>
half caObjH::dot(caObjH *source, int incx, int incy) {
  DEBUG_TRACE("Not implemented for half precision");
  return __float2half(0.f);
}
template <>
void caObjH::rot(caObjH *source, int incx, int incy, half, half) {
  DEBUG_TRACE("Not implemented for half precision");
}
template <>
void caObjH::gemv(char trans, half alpha, caObjH *matA, int lda, caObjH *vectx,
                  int incx, half beta, int incy) {
  int k = (((trans == 'N') || (trans == 'n')) ? matA->dims_data[2]
                                              : matA->dims_data[1]);

  carma_gemm(current_context->get_cublasHandle(), trans, 'N',
             this->dims_data[1], 1, k, alpha, matA->d_data, lda, vectx->d_data,
             vectx->dims_data[1], beta, this->d_data, this->dims_data[1]);
}
template <>
void caObjH::ger(half, caObjH *, int, caObjH *, int, int) {
  DEBUG_TRACE("Not implemented for half precision");
}
template <>
void caObjH::symv(char, half, caObjH *, int, caObjH *, int, half, int) {
  DEBUG_TRACE("Not implemented for half precision");
}
template void caObjH::gemm(char, char, half, caObjH *, int, caObjH *, int, half,
                           int);
template <>
void caObjH::symm(char, char, half, caObjH *, int, caObjH *, int, half, int) {
  DEBUG_TRACE("Not implemented for half precision");
}
template <>
void caObjH::syrk(char, char, half, caObjH *, int, half, int) {
  DEBUG_TRACE("Not implemented for half precision");
}
template <>
void caObjH::syrkx(char, char, half, caObjH *, int, caObjH *, int, half, int) {
  DEBUG_TRACE("Not implemented for half precision");
}
template void caObjH::geam(char, char, half, caObjH *, int, half, caObjH *, int,
                           int);
template void caObjH::dgmm(char, caObjH *, int, caObjH *, int, int);
template class carma_obj<half>;

#endif
