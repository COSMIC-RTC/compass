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

template caObjS::carma_obj(carma_context *current_context,
                           const long *dims_data);
template caObjD::carma_obj(carma_context *current_context,
                           const long *dims_data);
template caObjI::carma_obj(carma_context *current_context,
                           const long *dims_data);
template caObjUI::carma_obj(carma_context *current_context,
                            const long *dims_data);
template caObjS2::carma_obj(carma_context *current_context,
                            const long *dims_data);
template caObjD2::carma_obj(carma_context *current_context,
                            const long *dims_data);
template carma_obj<struct tuple_t<float> >::carma_obj(
    carma_context *current_context, const long *dims_data);

template <class T_data>
carma_obj<T_data>::carma_obj(carma_obj<T_data> *src) {
  /** \brief carma_obj creator.
   * \param src : carma_obj to copy
   */
  init(src->current_context, src->dims_data, src->d_data, false,
       src->get_nbStreams());
}
template caObjI::carma_obj(caObjI *src);
template caObjUI::carma_obj(caObjUI *src);
template caObjS::carma_obj(caObjS *src);
template caObjD::carma_obj(caObjD *src);
template caObjS2::carma_obj(caObjS2 *src);
template caObjD2::carma_obj(caObjD2 *src);

template <class T_data>
carma_obj<T_data>::carma_obj(carma_context *current_context,
                             carma_obj<T_data> *src) {
  /** \brief carma_obj creator.
   * \param src : carma_obj to copy
   */
  init(current_context, src->dims_data, src->d_data, false, 0);
}

template caObjI::carma_obj(carma_context *current_context, caObjI *src);
template caObjUI::carma_obj(carma_context *current_context, caObjUI *src);
template caObjS::carma_obj(carma_context *current_context, caObjS *src);
template caObjD::carma_obj(carma_context *current_context, caObjD *src);
template caObjS2::carma_obj(carma_context *current_context, caObjS2 *src);
template caObjD2::carma_obj(carma_context *current_context, caObjD2 *src);

template <class T_data>
carma_obj<T_data>::carma_obj(carma_context *current_context,
                             const long *dims_data, T_data *data) {
  /** \brief carma_obj creator.
   * \param dims_data : the array size Yorick format : [ndims,dims1,dims2,...]
   * \param data : the array
   */
  init(current_context, dims_data, data, true, 0);
}

template caObjS::carma_obj(carma_context *current_context,
                           const long *dims_data, float *data);
template caObjD::carma_obj(carma_context *current_context,
                           const long *dims_data, double *data);
template caObjI::carma_obj(carma_context *current_context,
                           const long *dims_data, int *data);
template caObjUI::carma_obj(carma_context *current_context,
                            const long *dims_data, unsigned int *data);
template caObjS2::carma_obj(carma_context *current_context,
                            const long *dims_data, float2 *data);
template caObjD2::carma_obj(carma_context *current_context,
                            const long *dims_data, double2 *data);

template <class T_data>
carma_obj<T_data>::carma_obj(carma_context *current_context,
                             const long *dims_data, int nb_streams) {
  /** \brief carma_obj creator.
   * \param dims_data : the array size Yorick format : [ndims,dims1,dims2,...]
   */
  init(current_context, dims_data, NULL, true, nb_streams);
}

template caObjS::carma_obj(carma_context *current_context,
                           const long *dims_data, int nb_streams);
template caObjD::carma_obj(carma_context *current_context,
                           const long *dims_data, int nb_streams);
template caObjI::carma_obj(carma_context *current_context,
                           const long *dims_data, int nb_streams);
template caObjUI::carma_obj(carma_context *current_context,
                            const long *dims_data, int nb_streams);
template caObjS2::carma_obj(carma_context *current_context,
                            const long *dims_data, int nb_streams);
template caObjD2::carma_obj(carma_context *current_context,
                            const long *dims_data, int nb_streams);

template <class T_data>
carma_obj<T_data>::carma_obj(carma_context *current_context,
                             const long *dims_data, T_data *data,
                             int nb_streams) {
  /** \brief carma_obj creator.
   * \param dims_data : the array size Yorick format : [ndims,dims1,dims2,...]
   * \param data : the array
   */
  init(current_context, dims_data, data, true, nb_streams);
}

template caObjS::carma_obj(carma_context *current_context,
                           const long *dims_data, float *data, int nb_streams);
template caObjD::carma_obj(carma_context *current_context,
                           const long *dims_data, double *data, int nb_streams);
template caObjI::carma_obj(carma_context *current_context,
                           const long *dims_data, int *data, int nb_streams);
template caObjUI::carma_obj(carma_context *current_context,
                            const long *dims_data, unsigned int *data,
                            int nb_streams);
template caObjS2::carma_obj(carma_context *current_context,
                            const long *dims_data, float2 *data,
                            int nb_streams);
template caObjD2::carma_obj(carma_context *current_context,
                            const long *dims_data, double2 *data,
                            int nb_streams);

template <class T_data>
void carma_obj<T_data>::init(carma_context *context, const long *dims_data,
                             T_data *data, bool fromHost, int nb_streams) {
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

  carmaSafeCall(cudaFree(this->d_data));
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

template caObjS::~carma_obj();
template caObjD::~carma_obj();
template caObjI::~carma_obj();
template caObjUI::~carma_obj();
template caObjD2::~carma_obj();
template caObjS2::~carma_obj();
template carma_obj<struct tuple_t<float> >::~carma_obj();

template <class T_data>
int carma_obj<T_data>::host2device(T_data *data) {
  /** \brief host2device data transfer.
   * \param data : input data
   *
   * this method fills d_input with the imput data
   */
  carmaSafeCall(cudaMemcpy(this->d_data, data, sizeof(T_data) * this->nb_elem,
                           cudaMemcpyHostToDevice));

  return EXIT_SUCCESS;
}

template int caObjI::host2device(int *data);
template int caObjUI::host2device(unsigned int *data);
template int caObjS::host2device(float *data);
template int caObjD::host2device(double *data);
template int caObjS2::host2device(float2 *data);
template int caObjD2::host2device(double2 *data);
template int carma_obj<struct tuple_t<float> >::host2device(
    struct tuple_t<float> *data);

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
int carma_obj<T_data>::device2host(T_data *data) {
  /** \brief device2host data transfer.
   * \param data : output data
   *
   * this method copies the values in d_output to the output array
   */
  carmaSafeCall(cudaMemcpy(data, this->d_data, sizeof(T_data) * this->nb_elem,
                           cudaMemcpyDeviceToHost));

  return EXIT_SUCCESS;
}

template int caObjS::device2host(float *data);
template int caObjD::device2host(double *data);
template int caObjI::device2host(int *data);
template int caObjUI::device2host(unsigned int *data);
template int caObjS2::device2host(float2 *data);
template int caObjD2::device2host(double2 *data);

template <class T_data>
int carma_obj<T_data>::host2deviceAsync(T_data *data, cudaStream_t stream) {
  /** \brief host2device data transfer.
   * \param data : input data
   *
   * this method fills d_input with the imput data
   */
  cudaMemcpyAsync(this->d_data, data, sizeof(T_data) * this->nb_elem,
                  cudaMemcpyHostToDevice, stream);

  return EXIT_SUCCESS;
}

template int caObjI::host2deviceAsync(int *data, cudaStream_t stream);
template int caObjUI::host2deviceAsync(unsigned int *data, cudaStream_t stream);
template int caObjS::host2deviceAsync(float *data, cudaStream_t stream);
template int caObjD::host2deviceAsync(double *data, cudaStream_t stream);
template int caObjS2::host2deviceAsync(float2 *data, cudaStream_t stream);
template int caObjD2::host2deviceAsync(double2 *data, cudaStream_t stream);

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

template int caObjS::device2hostAsync(float *data, cudaStream_t stream);
template int caObjD::device2hostAsync(double *data, cudaStream_t stream);
template int caObjI::device2hostAsync(int *data, cudaStream_t stream);
template int caObjUI::device2hostAsync(unsigned int *data, cudaStream_t stream);
template int caObjS2::device2hostAsync(float2 *data, cudaStream_t stream);
template int caObjD2::device2hostAsync(double2 *data, cudaStream_t stream);

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

template int caObjS::device2hostOpt(float *data);
template int caObjD::device2hostOpt(double *data);
template int caObjI::device2hostOpt(int *data);
template int caObjUI::device2hostOpt(unsigned int *data);
template int caObjS2::device2hostOpt(float2 *data);
template int caObjD2::device2hostOpt(double2 *data);

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
template int caObjS::copyInto(float *data, int nb_elem);
template int caObjD::copyInto(double *data, int nb_elem);
template int caObjC::copyInto(cuFloatComplex *data, int nb_elem);
template int caObjZ::copyInto(cuDoubleComplex *data, int nb_elem);
template int caObjI::copyInto(int *data, int nb_elem);
template int caObjUI::copyInto(unsigned int *data, int nb_elem);

template <class T_data>
int carma_obj<T_data>::copyFrom(T_data *data, int nb_elem) {
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
template int caObjS::copyFrom(float *data, int nb_elem);
template int caObjD::copyFrom(double *data, int nb_elem);
template int caObjC::copyFrom(cuFloatComplex *data, int nb_elem);
template int caObjZ::copyFrom(cuDoubleComplex *data, int nb_elem);
template int caObjI::copyFrom(int *data, int nb_elem);
template int caObjUI::copyFrom(unsigned int *data, int nb_elem);

template <class T_data>
T_data carma_obj<T_data>::sum() {
  int nBlocks;
  int nThreads;

  this->current_context->set_activeDevice(device, 1);
  sumGetNumBlocksAndThreads(this->nb_elem,
                            this->current_context->get_device(device), nBlocks,
                            nThreads);

  reduce<T_data>(this->nb_elem, nThreads, nBlocks, this->d_data, this->d_data);

  carmaCheckMsg("Kernel execution failed");

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
  T_data *h_odata = new T_data[nBlocks];
  cudaMemcpy(h_odata, this->d_data, sizeof(T_data), cudaMemcpyDeviceToHost);
  return h_odata[0];
}
template float caObjS::sum();
template double caObjD::sum();
/*add template cuFloatComplex
 caObjC::sum();*/

template <class T_data>
void carma_obj<T_data>::clip(T_data min, T_data max) {
  clip_array<T_data>(this->d_data, min, max, this->nb_elem,
                     this->current_context->get_device(device));
}
template void caObjS::clip(float min, float max);
template void caObjD::clip(double min, double max);

template <class T_data>
int carma_obj<T_data>::transpose(carma_obj<T_data> *source) {
  transposeCU(source->d_data, this->d_data, this->dims_data[1],
              this->dims_data[2]);
  return EXIT_SUCCESS;
}
template int caObjS::transpose(caObjS *data);
template int caObjD::transpose(caObjD *data);
template int caObjC::transpose(caObjC *data);
template int caObjZ::transpose(caObjZ *data);

/*
 __   __         _      _
 \ \ / /__  _ __(_) ___| | __ __      ___ __ __ _ _ __  _ __   ___ _ __ ___
 \ V / _ \| '__| |/ __| |/ / \ \ /\ / / '__/ _` | '_ \| '_ \ / _ \ '__/ __|
 | | (_) | |  | | (__|   <   \ V  V /| | | (_| | |_) | |_) |  __/ |  \__ \
  |_|\___/|_|  |_|\___|_|\_\   \_/\_/ |_|  \__,_| .__/| .__/ \___|_|  |___/
 |_|   |_|
 */

int _genericCUS(caObjS *dest, caObjS *source) {
  /*! \brief sum wrapper for Yorick (single).
   * \param handle : a pointer to the carma_obj object we want to sum
   *
   * this function returns the sum of the single-precision array
   */
  return launch_generic1d((float *)*source, (float *)*dest, source->getNbElem(),
                          dest->getContext()->get_device(dest->getDevice()));
}
int _genericCUD(caObjD *dest, caObjD *source) {
  /*! \brief sum wrapper for Yorick (double).
   * \param handle : a pointer to the carma_obj object we want to sum
   *
   * this function returns the sum of the double-precision array
   */
  return launch_generic1d((double *)*source, (double *)*dest,
                          source->getNbElem(),
                          dest->getContext()->get_device(dest->getDevice()));
}

int _generic2CUS(caObjS *dest, caObjS *source) {
  /*! \brief sum wrapper for Yorick (single).
   * \param handle : a pointer to the carma_obj object we want to sum
   *
   * this function returns the sum of the single-precision array
   */
  return launch_generic2d((float *)*source, (float *)*dest, source->getDims(1),
                          source->getDims(2));
}
int _generic2CUD(caObjD *dest, caObjD *source) {
  /*! \brief sum wrapper for Yorick (double).
   * \param handle : a pointer to the carma_obj object we want to sum
   *
   * this function returns the sum of the double-precision array
   */
  return launch_generic2d((double *)*source, (double *)*dest,
                          source->getDims(1), source->getDims(2));
}
