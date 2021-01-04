// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2019 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.
//  Distributed under GNU - LGPL
//
//  COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
//  General Public License as published by the Free Software Foundation, either version 3 of the License,
//  or any later version.
//
//  COMPASS: End-to-end AO simulation tool using GPU acceleration
//  The COMPASS platform was designed to meet the need of high-performance for the simulation of AO systems.
//
//  The final product includes a software package for simulating all the critical subcomponents of AO,
//  particularly in the context of the ELT and a real-time core based on several control approaches,
//  with performances consistent with its integration into an instrument. Taking advantage of the specific
//  hardware architecture of the GPU, the COMPASS tool allows to achieve adequate execution speeds to
//  conduct large simulation campaigns called to the ELT.
//
//  The COMPASS platform can be used to carry a wide variety of simulations to both testspecific components
//  of AO of the E-ELT (such as wavefront analysis device with a pyramid or elongated Laser star), and
//  various systems configurations such as multi-conjugate AO.
//
//  COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
//  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//  See the GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License along with COMPASS.
//  If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>.
// -----------------------------------------------------------------------------

//! \file      carma_obj.cpp
//! \ingroup   libcarma
//! \class     CarmaObj
//! \brief     this class provides wrappers to the generic carma object
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.0.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#include <carma_magma.h>
#include <carma_obj.h>
#include <cstdlib> /* required for randomize() and random() */

template <class T_data>
CarmaObj<T_data>::CarmaObj(CarmaContext *current_context,
                             const long *dims_data) {
  /** \brief CarmaObj creator.
   * \param current_context : the context in which the CarmaObj is created
   * \param dims_data : the array size Yorick format : [ndims,dims1,dims2,...]
   */
  init(current_context, dims_data, NULL, true, 0);
}

template <class T_data>
CarmaObj<T_data>::CarmaObj(const CarmaObj<T_data> *src) {
  /** \brief CarmaObj creator.
   * \param src : CarmaObj to copy
   */
  init(src->current_context, src->dims_data, src->d_data, false,
       src->get_nb_streams());
}

template <class T_data>
CarmaObj<T_data>::CarmaObj(CarmaContext *current_context,
                             const CarmaObj<T_data> *src) {
  /** \brief CarmaObj creator.
   * \param src : CarmaObj to copy
   */
  init(current_context, src->dims_data, src->d_data, false, 0);
}

template <class T_data>
CarmaObj<T_data>::CarmaObj(CarmaContext *current_context,
                             const long *dims_data, const T_data *data) {
  /** \brief CarmaObj creator.
   * \param dims_data : the array size Yorick format : [ndims,dims1,dims2,...]
   * \param data : the array
   */
  init(current_context, dims_data, data, true, 0);
}

template <class T_data>
CarmaObj<T_data>::CarmaObj(CarmaContext *current_context,
                             const long *dims_data, int nb_streams) {
  /** \brief CarmaObj creator.
   * \param dims_data : the array size Yorick format : [ndims,dims1,dims2,...]
   */
  init(current_context, dims_data, NULL, true, nb_streams);
}

template <class T_data>
CarmaObj<T_data>::CarmaObj(CarmaContext *current_context,
                             const long *dims_data, const T_data *data,
                             int nb_streams) {
  /** \brief CarmaObj creator.
   * \param dims_data : the array size Yorick format : [ndims,dims1,dims2,...]
   * \param data : the array
   */
  init(current_context, dims_data, data, true, nb_streams);
}

template <class T_data>
void CarmaObj<T_data>::init(CarmaContext *context, const long *dims_data,
                             const T_data *data, bool fromHost,
                             int nb_streams) {
  this->current_context = context;
  const long size_data = dims_data[0] + 1;
  this->dims_data = new long[size_data];
  memcpy(this->dims_data, dims_data, size_data * sizeof(long));

  this->nb_elem = dims_data[1];
  for (int i = 2; i < size_data; i++) this->nb_elem *= dims_data[i];

  carma_safe_call(
      cudaMalloc((void **)&(this->d_data), sizeof(T_data) * this->nb_elem));
  this->device = current_context->get_active_device();

  if (data == NULL)
    carma_safe_call(cudaMemset(this->d_data, 0, sizeof(T_data) * this->nb_elem));

  else if (fromHost)
    this->host2device(data);
  else
    this->copy_from(data, this->nb_elem);

  this->plan = 0;
  this->gen = 0;
  this->d_states = 0;
  this->values = 0;
  this->o_data = 0;
  carma_safe_call(cudaMalloc((void **)&(this->d_num_valid), sizeof(size_t)));

  this->streams = new CarmaStreams();
  this->add_stream(nb_streams);

#if DEBUG
  printf("CARMA Object created @ 0x%p on GPU%d\n", this,
         current_context->get_active_real_device());
#endif
}

/*
 template<class T_data>
 CarmaObj<T_data>& CarmaObj<T_data>::operator= (const CarmaObj<T_data>& obj){

 if(this->d_data!=0) carma_safe_call( cudaFree(this->d_data) );
 if(this->dims_data!=0) delete(this->dims_data);
 CarmaObj<T_data> *new_obj = new CarmaObj(CarmaContext *current_context,
 obj); return new_obj;
 }
 */

template <class T_data>
CarmaObj<T_data>::~CarmaObj() {
  /** \brief CarmaObj destructor.
   */

  int old_device = current_context->get_active_device();
  current_context->set_active_device(this->device, 1);

  dealloc();
  this->d_data = nullptr;

  delete[](this->dims_data);
  this->dims_data = 0;

  cudaFree(this->d_num_valid);
  this->d_num_valid = nullptr;

  delete this->streams;

  if (this->values != 0) carma_safe_call(cudaFree(this->values));

  // if (this->o_data!=0) carma_safe_call( cudaFree(this->o_data) );

  /* Destroy the CUFFT plan. */
  if (this->plan != 0) carmafft_safe_call(cufftDestroy(this->plan));

  if (this->gen != 0) this->destroy_prng_host();

  if (this->d_states != 0) this->destroy_prng();
    //  carma_safe_call(cudaFree(this->d_states));
    // this->d_states = 0;

#if DEBUG
  printf("CARMA Object deleted @ 0x%p on GPU%d\n", this,
         current_context->get_active_real_device());
#endif
  current_context->set_active_device(old_device, 1);
}

template <class T_data>
template <typename T_dest>
int CarmaObj<T_data>::host2device(const T_dest *data) {
  /** \brief host2device data transfer.
   * \param data : input data
   *
   * this method fills d_input with the imput data
   */

  carma_safe_call(cudaMemcpy(this->d_data, data, sizeof(T_dest) * this->nb_elem,
                           cudaMemcpyHostToDevice));

  return EXIT_SUCCESS;
}
template int CarmaObj<unsigned int>::template host2device<unsigned int>(
    const unsigned int *data);
template int CarmaObj<uint16_t>::template host2device<uint16_t>(
    const uint16_t *data);
template int CarmaObj<int>::template host2device<int>(const int *data);
template int CarmaObj<float>::template host2device<float>(const float *data);
template int CarmaObj<double>::template host2device<double>(
    const double *data);
template int CarmaObj<cuFloatComplex>::template host2device<cuFloatComplex>(
    const cuFloatComplex *data);
template int CarmaObj<cuDoubleComplex>::template host2device<cuDoubleComplex>(
    const cuDoubleComplex *data);

/*
 template<class T_data>
 T_data* CarmaObj<T_data>::get_data(){
 * \brief get_data data transfer.
 * \return data : pointer on the data
 *

 return d_data;
 }

 template float* CarmaObjS::get_data();
 template double* CarmaObjD::get_data();
 template unsigned int* CarmaObjUI::get_data();
 template float2* CarmaObjS2::get_data();
 template double2* CarmaObjD2::get_data();
 */

template <class T_data>
template <typename T_dest>
int CarmaObj<T_data>::device2host(T_dest *data) {
  /** \brief device2host data transfer.
   * \param data : output data
   *
   * this method copies the values in d_output to the output array
   */

  carma_safe_call(cudaMemcpy(data, this->d_data, sizeof(T_dest) * this->nb_elem,
                           cudaMemcpyDeviceToHost));

  return EXIT_SUCCESS;
}

template int CarmaObj<unsigned int>::template device2host<unsigned int>(
    unsigned int *data);
template int CarmaObj<int>::template device2host<int>(int *data);
template int CarmaObj<uint16_t>::template device2host<uint16_t>(
    uint16_t *data);
template int CarmaObj<float>::template device2host<float>(float *data);
template int CarmaObj<double>::template device2host<double>(double *data);
template int CarmaObj<cuFloatComplex>::template device2host<cuFloatComplex>(
    cuFloatComplex *data);
template int CarmaObj<cuDoubleComplex>::template device2host<cuDoubleComplex>(
    cuDoubleComplex *data);

#ifdef CAN_DO_HALF
template int CarmaObj<half>::template host2device<half>(const half *data);
template int CarmaObj<half>::template device2host<half>(half *data);
template <>
template <>
// std::enable_if_t<std::is_same<T_data, half>::value>
int CarmaObj<half>::host2device<float>(const float *data) {
  return copy_from_float_to_half(data, this->d_data, this->nb_elem,
                             this->current_context->get_device(this->device));
}
template <>
template <>
// std::enable_if_t<std::is_same<T_data, half>::value>
int CarmaObj<half>::device2host<float>(float *data) {
  return copy_from_half_to_float(this->d_data, data, this->nb_elem,
                             this->current_context->get_device(this->device));
}
#endif

template <class T_data>
int CarmaObj<T_data>::host2device_async(const T_data *data,
                                        cudaStream_t stream) {
  /** \brief host2device data transfer.
   * \param data : input data
   *
   * this method fills d_input with the imput data
   */

  carma_safe_call(cudaMemcpyAsync(this->d_data, data,
                                sizeof(T_data) * this->nb_elem,
                                cudaMemcpyHostToDevice, stream));

  return EXIT_SUCCESS;
}

/*
 template<class T_data>
 T_data* CarmaObj<T_data>::get_data(){
 * \brief get_data data transfer.
 * \return data : pointer on the data
 *

 return d_data;
 }

 template float* CarmaObjS::get_data();
 template double* CarmaObjD::get_data();
 template unsigned int* CarmaObjUI::get_data();
 template float2* CarmaObjS2::get_data();
 template double2* CarmaObjD2::get_data();
 */

template <class T_data>
int CarmaObj<T_data>::device2host_async(T_data *data, cudaStream_t stream) {
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
int CarmaObj<T_data>::device2host_opt(T_data *data) {
  /** \brief device2host data transfer.
   * \param data : output data
   *
   * this method copies the values in d_output to the output array
   */
  if (this->o_data == 0) return EXIT_FAILURE;

  carma_safe_call(cudaMemcpy(data, this->o_data, sizeof(T_data) * this->nb_elem,
                           cudaMemcpyDeviceToHost));

  return EXIT_SUCCESS;
}

template <class T_data>
int CarmaObj<T_data>::copy_into(T_data *data, int nb_elem) {
  /** \brief device2host data transfer.
   * \param data : output data
   *
   * this method copies the values in d_output to the output array
   */
  if (nb_elem > this->nb_elem) nb_elem = this->nb_elem;

  carma_safe_call(cudaMemcpy(data, this->d_data, sizeof(T_data) * nb_elem,
                           cudaMemcpyDeviceToDevice));

  return EXIT_SUCCESS;
}

template <class T_data>
int CarmaObj<T_data>::copy_from(const T_data *data, int nb_elem) {
  /** \brief device2host data transfer.
   * \param data : output data
   *
   * this method copies the values in d_output to the output array
   */
  if (nb_elem > this->nb_elem) nb_elem = this->nb_elem;

  carma_safe_call(cudaMemcpy(this->d_data, data, sizeof(T_data) * nb_elem,
                           cudaMemcpyDeviceToDevice));

  return EXIT_SUCCESS;
}

template <class T_data>
int CarmaObj<T_data>::copy_from_async(const T_data *data, int nb_elem, cudaStream_t stream) {
  /** \brief device2host data transfer.
   * \param data : output data
   *
   * this method copies the values in d_output to the output array
   */
  if (nb_elem > this->nb_elem) nb_elem = this->nb_elem;

  carma_safe_call(cudaMemcpyAsync(this->d_data, data, sizeof(T_data) * nb_elem,
                           cudaMemcpyDeviceToDevice, stream));

  return EXIT_SUCCESS;
}

#ifdef USE_OCTOPUS
template <class T_data>
int CarmaObj<T_data>::copy_into(ipc::Cacao<T_data> *cacaoInterface) {
  /** \brief device2host data transfer.
   * \param data : output data
   *
   * this method copies the values in d_output to the output array
   */
  this->copy_into((T_data *)cacaoInterface->outputPtr(), this->nb_elem);

  return EXIT_SUCCESS;
}

template <class T_data>
int CarmaObj<T_data>::copy_from(ipc::Cacao<T_data> *cacaoInterface) {
  /** \brief device2host data transfer.
   * \param data : output data
   *
   * this method copies the values in d_output to the output array
   */
  this->copy_from((T_data *)cacaoInterface->outputPtr(), this->nb_elem);

  return EXIT_SUCCESS;
}

#endif
template <class T_data>
T_data CarmaObj<T_data>::sum() {
  return reduce<T_data>(this->d_data, this->nb_elem);
  // int nb_blocks;
  // int nb_threads;

  // this->current_context->set_active_device(device, 1);
  // sum_get_num_blocks_and_threads(this->nb_elem,
  //                           this->current_context->get_device(device),
  //                           nb_blocks, nb_threads);

  // reduce<T_data>(this->nb_elem, nb_threads, nb_blocks, this->d_data,
  // this->d_data);

  // carma_check_msg("Kernel execution failed");

  // sum partial block sums on GPU
  /*
  int s = nb_blocks;
  while (s > 1) {
    int threads = 0, blocks = 0;
    this->current_context->set_active_device(device,1);
    sum_get_num_blocks_and_threads(s, this->current_context->get_device(device),
  blocks, threads);

    reduce<T_data>(s, threads, blocks, this->d_data, this->d_data);

    s = (s + (threads * 2 - 1)) / (threads * 2);

  }
  */
  // T_data *h_odata = new T_data[nb_blocks];
  // cudaMemcpy(h_odata, this->d_data, sizeof(T_data), cudaMemcpyDeviceToHost);
  // return h_odata[0];
}

template <class T_data>
void CarmaObj<T_data>::init_reduceCub() {
  init_reduceCubCU(this->cub_data, this->cub_data_size, this->d_data,
                   this->o_data, this->nb_elem);
}

template <class T_data>
void CarmaObj<T_data>::reduceCub(cudaStream_t stream) {
  reduceCubCU(this->cub_data, this->cub_data_size, this->d_data, this->o_data,
              this->nb_elem, stream);
}

template <class T_data>
void CarmaObj<T_data>::clip(T_data min, T_data max, cudaStream_t stream) {
  clip_array<T_data>(this->d_data, min, max, this->nb_elem,
                     this->current_context->get_device(device), stream);
}

template <class T_data>
int CarmaObj<T_data>::transpose(CarmaObj<T_data> *source) {
  transposeCU(source->d_data, this->d_data, this->dims_data[1],
              this->dims_data[2]);
  return EXIT_SUCCESS;
}

template <class T>
int CarmaObj<T>::init_prng(long seed) {
  /*
    cudaDeviceProp deviceProperties =
    current_context->get_device(device)->get_properties(); const int maxThreads
    = deviceProperties.maxThreadsPerBlock; const int maxBlockDim =
    (deviceProperties.maxThreadsDim)[0]; int genPerBlock = std::min(maxThreads,
    maxBlockDim) / 2; int blockCount = deviceProperties.multiProcessorCount * 2;
  */
  int genPerBlock, blockCount;
  get_num_blocks_and_threads(current_context->get_device(device), nb_elem,
                         blockCount, genPerBlock);

  // Allocate memory for RNG states
  if (this->d_states == NULL)
    carma_safe_call(cudaMalloc((void **)&(this->d_states),
                             blockCount * genPerBlock * sizeof(curandState)));
  carma_safe_call(cudaMemset(this->d_states, 0,
                           blockCount * genPerBlock * sizeof(curandState)));

  this->nb_threads = genPerBlock;
  this->nb_blocks = blockCount;
  // randomize();
  std::vector<int> aseed(genPerBlock * blockCount);
  for (int cc = 0; cc <= genPerBlock * blockCount; cc++)
    aseed[cc] = seed + cc;  // random();

  int *seeds;
  carma_safe_call(
      cudaMalloc((void **)&seeds, genPerBlock * blockCount * sizeof(int)));
  carma_safe_call(cudaMemcpy(seeds, aseed.data(),
                           genPerBlock * blockCount * sizeof(int),
                           cudaMemcpyHostToDevice));

  // cerr << genPerBlock << " | " << blockCount << endl;
  carma_prng_init(seeds, genPerBlock, blockCount, this->d_states);

  carma_safe_call(cudaFree(seeds));

  return EXIT_SUCCESS;
}

#ifdef CAN_DO_HALF
template <>
int CarmaObjH::init_prng(long seed) {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_SUCCESS;
}
#endif

template <class T>
int CarmaObj<T>::init_prng() {
  return this->init_prng(1234);
}

#ifdef CAN_DO_HALF
template <>
int CarmaObjH::init_prng() {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_SUCCESS;
}
#endif

template <class T>
int CarmaObj<T>::destroy_prng() {
  carma_safe_call(cudaFree(this->d_states));
  return EXIT_SUCCESS;
}

#ifdef CAN_DO_HALF
template <>
int CarmaObjH::destroy_prng() {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_SUCCESS;
}
#endif

template <class T>
int CarmaObj<T>::prng(T *output, char gtype, float alpha, float beta) {
  carma_prng_cu(output, this->nb_threads, this->nb_blocks, this->d_states, gtype,
                this->nb_elem, alpha, beta);

  return EXIT_SUCCESS;
}

#ifdef CAN_DO_HALF
template <>
int CarmaObjH::prng(half *output, char gtype, float alpha, float beta) {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_SUCCESS;
}
#endif

template <class T>
int CarmaObj<T>::prng_montagn(float init_montagn) {
  carma_curand_montagn(this->d_states, this->d_data, this->nb_elem,
                       this->current_context->get_device(this->device));

  return EXIT_SUCCESS;
}

#ifdef CAN_DO_HALF
template <>
int CarmaObjH::prng_montagn(float init_montagn) {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_SUCCESS;
}
#endif

template <class T>
int CarmaObj<T>::prng(T *output, char gtype, float alpha) {
  carma_prng_cu(output, this->nb_threads, this->nb_blocks, this->d_states, gtype,
                this->nb_elem, alpha, 0.0f);

  return EXIT_SUCCESS;
}

#ifdef CAN_DO_HALF
template <>
int CarmaObjH::prng(half *output, char gtype, float alpha) {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_SUCCESS;
}
#endif

template <class T>
int CarmaObj<T>::prng(char gtype) {
  return prng(this->d_data, gtype, 1.0f, 0.0f);
}

#ifdef CAN_DO_HALF
template <>
int CarmaObjH::prng(char gtype) {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_SUCCESS;
}
#endif

template <class T>
int CarmaObj<T>::prng(char gtype, float alpha) {
  return prng(this->d_data, gtype, alpha, 0.0f);
}

#ifdef CAN_DO_HALF
template <>
int CarmaObjH::prng(char gtype, float alpha) {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_SUCCESS;
}
#endif

template <class T>
int CarmaObj<T>::prng(char gtype, float alpha, float beta) {
  return prng(this->d_data, gtype, alpha, beta);
}

#ifdef CAN_DO_HALF
template <>
int CarmaObjH::prng(char gtype, float alpha, float beta) {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_SUCCESS;
}
#endif

template <class T>
int CarmaObj<T>::init_prng_host(int seed) {
  if (this->gen != NULL) curandDestroyGenerator(this->gen);

  curandCreateGenerator(&(this->gen), CURAND_RNG_PSEUDO_XORWOW);
  // CURAND_RNG_PSEUDO_MTGP32
  // CURAND_RNG_PSEUDO_XORWOW
  curandSetPseudoRandomGeneratorSeed(this->gen, seed);

  return EXIT_SUCCESS;
}

#ifdef CAN_DO_HALF
template <>
int CarmaObjH::init_prng_host(int seed) {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_SUCCESS;
}
#endif

template <class T>
int CarmaObj<T>::prng_host(char gtype) {
  DEBUG_TRACE("Not implemented");
  return EXIT_FAILURE;
}

template <>
int CarmaObjS::prng_host(char gtype) {
  if (gtype == 'U')
    curandGenerateUniform(this->gen, this->d_data, this->nb_elem);
  if (gtype == 'N')
    curandGenerateNormal(this->gen, this->d_data, this->nb_elem, 0.0f, 1.0f);
  return EXIT_SUCCESS;
}

template <>
int CarmaObjD::prng_host(char gtype) {
  if (gtype == 'U')
    curandGenerateUniformDouble(this->gen, this->d_data, this->nb_elem);
  if (gtype == 'N')
    curandGenerateNormalDouble(this->gen, this->d_data, this->nb_elem, 0.0,
                               1.0);
  return EXIT_SUCCESS;
}

template <class T>
int CarmaObj<T>::prng_host(char gtype, T stddev) {
  DEBUG_TRACE("Not implemented");
  return EXIT_FAILURE;
}

template <>
int CarmaObjS::prng_host(char gtype, float stddev) {
  if (gtype == 'U')
    curandGenerateUniform(this->gen, this->d_data, this->nb_elem);
  if (gtype == 'N')
    curandGenerateNormal(this->gen, this->d_data, this->nb_elem, 0.0f, stddev);
  return EXIT_SUCCESS;
}
template <>
int CarmaObjD::prng_host(char gtype, double stddev) {
  if (gtype == 'U')
    curandGenerateUniformDouble(this->gen, this->d_data, this->nb_elem);
  if (gtype == 'N')
    curandGenerateNormalDouble(this->gen, this->d_data, this->nb_elem, 0.0,
                               stddev);
  return EXIT_SUCCESS;
}

template <class T>
int CarmaObj<T>::prng_host(char gtype, T stddev, T alpha) {
  DEBUG_TRACE("Not implemented");
  return EXIT_FAILURE;
}
template <>
int CarmaObjS::prng_host(char gtype, float stddev, float alpha) {
  float *tmp_data;

  carma_safe_call(cudaMalloc((void **)&tmp_data, this->nb_elem * sizeof(float)));

  if (gtype == 'U') curandGenerateUniform(this->gen, tmp_data, this->nb_elem);
  if (gtype == 'N')
    curandGenerateNormal(this->gen, tmp_data, this->nb_elem, 0.0f, stddev);

  carma_axpy(current_context->get_cublas_handle(), this->nb_elem, alpha,
             tmp_data, 1, this->d_data, 1);

  carma_safe_call(cudaFree(tmp_data));

  return EXIT_SUCCESS;
}

template <>
int CarmaObjD::prng_host(char gtype, double stddev, double alpha) {
  double *tmp_data;

  carma_safe_call(cudaMalloc((void **)&tmp_data, this->nb_elem * sizeof(double)));

  if (gtype == 'U')
    curandGenerateUniformDouble(this->gen, tmp_data, this->nb_elem);
  if (gtype == 'N')
    curandGenerateNormalDouble(this->gen, tmp_data, this->nb_elem, 0.0f,
                               stddev);

  carma_axpy(current_context->get_cublas_handle(), this->nb_elem, alpha,
             tmp_data, 1, this->d_data, 1);

  carma_safe_call(cudaFree(tmp_data));

  return EXIT_SUCCESS;
}

template <class T>
int CarmaObj<T>::destroy_prng_host() {
  curandDestroyGenerator(this->gen);
  return EXIT_SUCCESS;
}

#ifdef CAN_DO_HALF
template <>
int CarmaObjH::destroy_prng_host() {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_SUCCESS;
}
#endif

template <class T_data>
int CarmaObj<T_data>::host2device_vect(const T_data *data, int incx, int incy) {
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
int CarmaObj<T_data>::device2host_vect(T_data *data, int incx, int incy) {
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
int CarmaObj<T_data>::host2device_mat(const T_data *data, int lda, int ldb) {
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
int CarmaObj<T_data>::device2host_mat(T_data *data, int lda, int ldb) {
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
int CarmaObj<T_data>::aimax(int incx) {
  /** \brief aimax method
   * \param incx : increment in x
   *
   * this method finds the smallest index of the maximum magnitude element in
   * vect
   */
  return carma_where_amax(current_context->get_cublas_handle(), this->nb_elem,
                          this->d_data, incx);
}

template <class T_data>
int CarmaObj<T_data>::aimin(int incx) {
  /** \brief aimin method
   * \param incx : increment in x
   *
   * this method finds the smallest index of the minimum magnitude element in
   * vect
   */
  return carma_where_amin(current_context->get_cublas_handle(), this->nb_elem,
                          this->d_data, incx);
};
template <class T_data>
T_data CarmaObj<T_data>::asum(int incx) {
  /** \brief asum method
   * \param incx : increment in x
   *
   * this method computes the sum of the absolute values of the elements
   */
  return carma_getasum(current_context->get_cublas_handle(), this->nb_elem,
                       this->d_data, incx);
  ;
};
template <class T_data>
T_data CarmaObj<T_data>::nrm2(int incx) {
  /** \brief getNrm2 method
   * \param n    : vect size
   * \param vect : vector to sum (x)
   * \param incx : increment in x
   *
   * this method computes the Euclidean norm of vect
   */
  return carma_nrm2(current_context->get_cublas_handle(), this->nb_elem,
                    this->d_data, incx);
}
template <class T_data>
T_data CarmaObj<T_data>::dot(CarmaObj<T_data> *source, int incx, int incy) {
  /** \brief ccomputeDot method
   * \param n    : vect size
   * \param vectx : first vector (x)
   * \param incx : increment in x
   * \param vecty : second vector (y)
   * \param incy : increment in y
   *
   * this method computes the dot product of two vectors
   */
  return carma_dot(current_context->get_cublas_handle(), this->nb_elem,
                   source->d_data, incx, this->d_data, incy);
}
template <class T_data>
void CarmaObj<T_data>::scale(T_data alpha, int incx) {
  /** \brief vectScale method
   * \param n    : vect size
   * \param alpha : scale factor
   * \param vect : vector to scale (x)
   * \param incx : increment in x
   *
   * this method replaces vector x with alpha * x

   */
  carma_scal(current_context->get_cublas_handle(), this->nb_elem, alpha,
             this->d_data, incx);
}
template <class T_data>
void CarmaObj<T_data>::swap(CarmaObj<T_data> *source, int incx, int incy) {
  /** \brief vectSwap method
   * \param n    : vect size
   * \param vectx : first vector to swap (x)
   * \param incx : increment in x
   * \param vecty : second vector to swap (y)
   * \param incy : increment in y
   *
   * this method interchanges vector x with vector y
   */

  carma_swap(current_context->get_cublas_handle(), this->nb_elem, source->d_data,
             incx, this->d_data, incy);
}

template <class T_data>
void CarmaObj<T_data>::copy(CarmaObj<T_data> *source, int incx, int incy) {
  /** \brief vectCopy method
   * \param n    : vect size
   * \param vectx : vector to copy (x)
   * \param incx : increment in x
   * \param vecty : vector to fill (y)
   * \param incy : increment in y
   *
   * this method interchanges vector x with vector y
   */

  carma_copy(current_context->get_cublas_handle(), this->nb_elem, source->d_data,
             incx, this->d_data, incy);
};
template <class T_data>
void CarmaObj<T_data>::axpy(T_data alpha, CarmaObj<T_data> *source, int incx,
                             int incy, int offset) {
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
  carma_axpy(current_context->get_cublas_handle(), this->nb_elem, alpha,
             source->get_data_at(offset), incx, this->d_data, incy);
}
template <class T_data>
void CarmaObj<T_data>::rot(CarmaObj<T_data> *source, int incx, int incy,
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
  carma_rot(current_context->get_cublas_handle(), this->nb_elem, source->d_data,
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
void CarmaObj<T_data>::gemv(char trans, T_data alpha, CarmaObj<T_data> *matA,
                             int lda, CarmaObj<T_data> *vectx, int incx,
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
  carma_gemv(current_context->get_cublas_handle(), trans, matA->dims_data[1],
             matA->dims_data[2], alpha, matA->d_data, lda, vectx->d_data, incx,
             beta, this->d_data, incy);
}
template <class T_data>
void CarmaObj<T_data>::ger(T_data alpha, CarmaObj<T_data> *vectx, int incx,
                            CarmaObj<T_data> *vecty, int incy, int lda) {
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
  carma_ger(current_context->get_cublas_handle(), this->dims_data[1],
            this->dims_data[2], alpha, vectx->d_data, incx, vecty->d_data, incy,
            this->d_data, lda);
}

template <class T_data>
void CarmaObj<T_data>::symv(char uplo, T_data alpha, CarmaObj<T_data> *matA,
                             int lda, CarmaObj<T_data> *vectx, int incx,
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
  carma_symv(current_context->get_cublas_handle(), uplo, matA->dims_data[1],
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
void CarmaObj<T_data>::gemm(char transa, char transb, T_data alpha,
                             CarmaObj<T_data> *matA, int lda,
                             CarmaObj<T_data> *matB, int ldb, T_data beta,
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
  carma_gemm(current_context->get_cublas_handle(), transa, transb,
             this->dims_data[1], this->dims_data[2], k, alpha, matA->d_data,
             lda, matB->d_data, ldb, beta, this->d_data, ldc);
}
template <class T_data>
void CarmaObj<T_data>::symm(char side, char uplo, T_data alpha,
                             CarmaObj<T_data> *matA, int lda,
                             CarmaObj<T_data> *matB, int ldb, T_data beta,
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
  carma_symm(current_context->get_cublas_handle(), side, uplo,
             this->dims_data[1], this->dims_data[2], alpha, matA->d_data, lda,
             matB->d_data, ldb, beta, this->d_data, ldc);
}
template <class T_data>
void CarmaObj<T_data>::syrk(char uplo, char transa, T_data alpha,
                             CarmaObj<T_data> *matA, int lda, T_data beta,
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
   carma_syrk(current_context->get_cublas_handle(), uplo, transa,
   this->dims_data[1], this->dims_data[2], alpha, matA->d_data, lda, beta,
   this->d_data, ldc);
   */
  carma_syrk(current_context->get_cublas_handle(), uplo, transa,
             matA->dims_data[1], matA->dims_data[2], alpha, matA->d_data, lda,
             beta, this->d_data, ldc);
  // "this" refer to matrix "C": this->dims_data[1]=this->dims_data[2]=n
}
template <class T_data>
void CarmaObj<T_data>::syrkx(char uplo, char transa, T_data alpha,
                              CarmaObj<T_data> *matA, int lda,
                              CarmaObj<T_data> *matB, int ldb, T_data beta,
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
  /*carma_syrkx(current_context->get_cublas_handle(), uplo, transa,
   this->dims_data[1], this->dims_data[2], alpha, matA->d_data, lda,
   matB->d_data, ldb, beta, this->d_data, ldc);*/
  carma_syrkx(current_context->get_cublas_handle(), uplo, transa,
              matA->dims_data[1], matA->dims_data[2], alpha, matA->d_data, lda,
              matB->d_data, ldb, beta, this->d_data, ldc);
}
template <class T_data>
void CarmaObj<T_data>::geam(char transa, char transb, T_data alpha,
                             CarmaObj<T_data> *matA, int lda, T_data beta,
                             CarmaObj<T_data> *matB, int ldb, int ldc) {
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

  carma_geam(current_context->get_cublas_handle(), transa, transb,
             this->dims_data[1], this->dims_data[2], alpha, matA->d_data, lda,
             beta, matB->d_data, ldb, this->d_data, ldc);
}
template <class T_data>
void CarmaObj<T_data>::dgmm(char side, CarmaObj<T_data> *matA, int lda,
                             CarmaObj<T_data> *vectx, int incx, int ldc) {
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

  carma_dgmm(current_context->get_cublas_handle(), side, this->dims_data[1],
             this->dims_data[2], matA->d_data, lda, vectx->d_data, incx,
             this->d_data, ldc);
}

template class CarmaObj<int>;
template class CarmaObj<unsigned int>;
template class CarmaObj<uint16_t>;
template class CarmaObj<float>;
template class CarmaObj<double>;
template class CarmaObj<cuFloatComplex>;
template class CarmaObj<cuDoubleComplex>;
// template class CarmaObj<struct tuple_t<float>>;

#ifdef CAN_DO_HALF
template <>
int CarmaObjH::aimax(int incx) {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_FAILURE;
}
template <>
int CarmaObjH::aimin(int incx) {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_FAILURE;
}
template <>
half CarmaObjH::asum(int incx) {
  DEBUG_TRACE("Not implemented for half precision");
  return __float2half(0.f);
}
template <>
half CarmaObjH::nrm2(int incx) {
  DEBUG_TRACE("Not implemented for half precision");
  return __float2half(0.f);
}
template <>
void CarmaObjH::scale(half alpha, int incx) {
  DEBUG_TRACE("Not implemented for half precision");
}
template <>
void CarmaObjH::copy(CarmaObjH *, int incx, int incy) {
  DEBUG_TRACE("Not implemented for half precision");
}
template <>
void CarmaObjH::swap(CarmaObjH *, int incx, int incy) {
  DEBUG_TRACE("Not implemented for half precision");
}
template <>
void CarmaObjH::axpy(half alpha, CarmaObjH *source, int incx, int incy, int offset) {
  custom_half_axpy(alpha, source->get_data_at(offset), incx, incy,
                   this->get_nb_elements(), this->d_data,
                   current_context->get_device(this->device));
}
template <>
half CarmaObjH::dot(CarmaObjH *source, int incx, int incy) {
  DEBUG_TRACE("Not implemented for half precision");
  return __float2half(0.f);
}
template <>
void CarmaObjH::rot(CarmaObjH *source, int incx, int incy, half, half) {
  DEBUG_TRACE("Not implemented for half precision");
}
template <>
void CarmaObjH::gemv(char trans, half alpha, CarmaObjH *matA, int lda, CarmaObjH *vectx,
                  int incx, half beta, int incy) {
  int k = (((trans == 'N') || (trans == 'n')) ? matA->dims_data[2]
                                              : matA->dims_data[1]);

  carma_gemm(current_context->get_cublas_handle(), trans, 'N',
             this->dims_data[1], 1, k, alpha, matA->d_data, lda, vectx->d_data,
             vectx->dims_data[1], beta, this->d_data, this->dims_data[1]);
}
template <>
void CarmaObjH::ger(half, CarmaObjH *, int, CarmaObjH *, int, int) {
  DEBUG_TRACE("Not implemented for half precision");
}
template <>
void CarmaObjH::symv(char, half, CarmaObjH *, int, CarmaObjH *, int, half, int) {
  DEBUG_TRACE("Not implemented for half precision");
}
template void CarmaObjH::gemm(char, char, half, CarmaObjH *, int, CarmaObjH *, int, half,
                           int);
template <>
void CarmaObjH::symm(char, char, half, CarmaObjH *, int, CarmaObjH *, int, half, int) {
  DEBUG_TRACE("Not implemented for half precision");
}
template <>
void CarmaObjH::syrk(char, char, half, CarmaObjH *, int, half, int) {
  DEBUG_TRACE("Not implemented for half precision");
}
template <>
void CarmaObjH::syrkx(char, char, half, CarmaObjH *, int, CarmaObjH *, int, half, int) {
  DEBUG_TRACE("Not implemented for half precision");
}
template void CarmaObjH::geam(char, char, half, CarmaObjH *, int, half, CarmaObjH *, int,
                           int);
template void CarmaObjH::dgmm(char, CarmaObjH *, int, CarmaObjH *, int, int);
template class CarmaObj<half>;

#endif
