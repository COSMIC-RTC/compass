#include <carma_obj.h>

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

template class carma_obj<int>;
template class carma_obj<unsigned int>;
template class carma_obj<float>;
template class carma_obj<double>;
template class carma_obj<cuFloatComplex>;
template class carma_obj<cuDoubleComplex>;
template class carma_obj<struct tuple_t<float>>;
#ifdef CAN_DO_HALF
template class carma_obj<half>;
#endif
