// This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
//
// COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU
// Lesser General Public License as published by the Free Software Foundation, either version 3 of
// the License, or any later version.
//
// COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
// even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License along with COMPASS. If
// not, see <https://www.gnu.org/licenses/>
//
//  Copyright (C) 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>

//! \file      carma_host_obj.cpp
//! \ingroup   libcarma
//! \class     CarmaHostObj
//! \brief     this class provides wrappers to the generic carma host object
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24


#include <carma_host_obj.hpp>
#include <carma_obj.hpp>

#include <cuda_runtime.h>
#include <sys/mman.h>  // for mmap() / munmap()

template <class T_data>
CarmaHostObj<T_data>::CarmaHostObj(const int64_t *dims_data) {
  /** \brief CarmaHostObj creator.
   * \param dims_data : the array size Yorick format : [ndims,dims1,dims2,...]
   */
  init(dims_data, 0L, MA_MALLOC, 0);
}
template CarmaHostObj<float>::CarmaHostObj(const int64_t *dims_data);
template CarmaHostObj<double>::CarmaHostObj(const int64_t *dims_data);
template CarmaHostObj<int32_t>::CarmaHostObj(const int64_t *dims_data);
template CarmaHostObj<uint32_t>::CarmaHostObj(const int64_t *dims_data);
template CarmaHostObj<cuFloatComplex>::CarmaHostObj(const int64_t *dims_data);
template CarmaHostObj<cuDoubleComplex>::CarmaHostObj(const int64_t *dims_data);


template <class T_data>
CarmaHostObj<T_data>::CarmaHostObj() {
  /** \brief CarmaHostObj creator.
   * \param dims_data : the array size Yorick format : [ndims,dims1,dims2,...]
   */

  int64_t dims_data[2] = {1,0};
  init(dims_data, 0L, MA_MALLOC, 0);
}

template <class T_data>
CarmaHostObj<T_data>::CarmaHostObj(const std::vector<int64_t> &dims):
CarmaHostObj(dims.data()){}

template <class T_data>
CarmaHostObj<T_data>::CarmaHostObj(const int64_t *dims_data,
                                       MemAlloc malloc_type) {
  /** \brief CarmaHostObj creator.
   * \param dims_data : the array size Yorick format : [ndims,dims1,dims2,...]
   * \param malloc_type : type of memory allacation : MA_MALLOC, MA_PAGELOCK,
   * MA_ZERO, MA_PORT, MA_WC
   */
  init(dims_data, 0L, malloc_type, 0);
}
template CarmaHostObj<float>::CarmaHostObj(const int64_t *dims_data,
                                               MemAlloc malloc_type);
template CarmaHostObj<double>::CarmaHostObj(const int64_t *dims_data,
                                                MemAlloc malloc_type);
template CarmaHostObj<int32_t>::CarmaHostObj(const int64_t *dims_data,
                                             MemAlloc malloc_type);
template CarmaHostObj<uint32_t>::CarmaHostObj(const int64_t *dims_data,
                                                      MemAlloc malloc_type);
template CarmaHostObj<cuFloatComplex>::CarmaHostObj(const int64_t *dims_data,
                                                        MemAlloc malloc_type);
template CarmaHostObj<cuDoubleComplex>::CarmaHostObj(const int64_t *dims_data,
                                                         MemAlloc malloc_type);

template <class T_data>
CarmaHostObj<T_data>::CarmaHostObj(const CarmaHostObj<T_data> *src) {
  /** \brief CarmaHostObj creator.
   * \param src : CarmaHostObj to copy
   */
  init(src->dims_data, src->h_data, src->malloc_type, 0);
}
template CarmaHostObj<float>::CarmaHostObj(
    const CarmaHostObj<float> *src);
template CarmaHostObj<double>::CarmaHostObj(
    const CarmaHostObj<double> *src);
template CarmaHostObj<int32_t>::CarmaHostObj(const CarmaHostObj<int32_t> *src);
template CarmaHostObj<uint32_t>::CarmaHostObj(
    const CarmaHostObj<uint32_t> *src);
template CarmaHostObj<cuFloatComplex>::CarmaHostObj(
    const CarmaHostObj<cuFloatComplex> *src);
template CarmaHostObj<cuDoubleComplex>::CarmaHostObj(
    const CarmaHostObj<cuDoubleComplex> *src);

template <class T_data>
CarmaHostObj<T_data>::CarmaHostObj(const CarmaHostObj<T_data> *src,
                                       MemAlloc malloc_type) {
  /** \brief CarmaHostObj creator.
   * \param src : CarmaHostObj to copy
   * \param malloc_type : type of memory allacation : MA_MALLOC, MA_PAGELOCK,
   * MA_ZERO, MA_PORT, MA_WC
   */
  init(src->dims_data, src->h_data, malloc_type, 0);
}
template CarmaHostObj<float>::CarmaHostObj(const CarmaHostObj<float> *src,
                                               MemAlloc malloc_type);
template CarmaHostObj<double>::CarmaHostObj(
    const CarmaHostObj<double> *src, MemAlloc malloc_type);
template CarmaHostObj<int32_t>::CarmaHostObj(const CarmaHostObj<int32_t> *src,
                                             MemAlloc malloc_type);
template CarmaHostObj<uint32_t>::CarmaHostObj(
    const CarmaHostObj<uint32_t> *src, MemAlloc malloc_type);
template CarmaHostObj<cuFloatComplex>::CarmaHostObj(
    const CarmaHostObj<cuFloatComplex> *src, MemAlloc malloc_type);
template CarmaHostObj<cuDoubleComplex>::CarmaHostObj(
    const CarmaHostObj<cuDoubleComplex> *src, MemAlloc malloc_type);

template <class T_data>
CarmaHostObj<T_data>::CarmaHostObj(const int64_t *dims_data,
                                       const T_data *data) {
  /** \brief CarmaHostObj creator.
   * \param dims_data : the array size Yorick format : [ndims,dims1,dims2,...]
   * \param data : the array
   */
  init(dims_data, data, MA_MALLOC, 0);
}
template CarmaHostObj<float>::CarmaHostObj(const int64_t *dims_data,
                                               const float *data);
template CarmaHostObj<double>::CarmaHostObj(const int64_t *dims_data,
                                                const double *data);
template CarmaHostObj<int32_t>::CarmaHostObj(const int64_t *dims_data,
                                             const int32_t *data);
template CarmaHostObj<uint32_t>::CarmaHostObj(const int64_t *dims_data,
                                                      const uint32_t *data);
template CarmaHostObj<cuFloatComplex>::CarmaHostObj(
    const int64_t *dims_data, const cuFloatComplex *data);
template CarmaHostObj<cuDoubleComplex>::CarmaHostObj(
    const int64_t *dims_data, const cuDoubleComplex *data);

template <class T_data>
CarmaHostObj<T_data>::CarmaHostObj(const int64_t *dims_data,
                                       const T_data *data,
                                       MemAlloc malloc_type) {
  /** \brief CarmaHostObj creator.
   * \param dims_data : the array size Yorick format : [ndims,dims1,dims2,...]
   * \param data : the array
   * \param malloc_type : type of memory allacation : MA_MALLOC, MA_PAGELOCK,
   * MA_ZERO, MA_PORT, MA_WC
   */
  init(dims_data, data, malloc_type, 0);
}
template CarmaHostObj<float>::CarmaHostObj(const int64_t *dims_data,
                                               const float *data,
                                               MemAlloc malloc_type);
template CarmaHostObj<double>::CarmaHostObj(const int64_t *dims_data,
                                                const double *data,
                                                MemAlloc malloc_type);
template CarmaHostObj<int32_t>::CarmaHostObj(const int64_t *dims_data,
                                             const int32_t *data,
                                             MemAlloc malloc_type);
template CarmaHostObj<uint32_t>::CarmaHostObj(const int64_t *dims_data,
                                                      const uint32_t *data,
                                                      MemAlloc malloc_type);
template CarmaHostObj<cuFloatComplex>::CarmaHostObj(
    const int64_t *dims_data, const cuFloatComplex *data, MemAlloc malloc_type);
template CarmaHostObj<cuDoubleComplex>::CarmaHostObj(
    const int64_t *dims_data, const cuDoubleComplex *data, MemAlloc malloc_type);

template <class T_data>
CarmaHostObj<T_data>::CarmaHostObj(const int64_t *dims_data, int32_t nb_streams) {
  /** \brief CarmaHostObj creator.
   * \param dims_data : the array size Yorick format : [ndims,dims1,dims2,...]
   */
  init(dims_data, 0L, MA_MALLOC, nb_streams);
}
template CarmaHostObj<float>::CarmaHostObj(const int64_t *dims_data,
                                               int32_t nb_streams);
template CarmaHostObj<double>::CarmaHostObj(const int64_t *dims_data,
                                                int32_t nb_streams);
template CarmaHostObj<int32_t>::CarmaHostObj(const int64_t *dims_data,
                                             int32_t nb_streams);
template CarmaHostObj<uint32_t>::CarmaHostObj(const int64_t *dims_data,
                                                      int32_t nb_streams);
template CarmaHostObj<cuFloatComplex>::CarmaHostObj(const int64_t *dims_data,
                                                        int32_t nb_streams);
template CarmaHostObj<cuDoubleComplex>::CarmaHostObj(const int64_t *dims_data,
                                                         int32_t nb_streams);

template <class T_data>
CarmaHostObj<T_data>::CarmaHostObj(const int64_t *dims_data,
                                       MemAlloc malloc_type, int32_t nb_streams) {
  /** \brief CarmaHostObj creator.
   * \param dims_data : the array size Yorick format : [ndims,dims1,dims2,...]
   * \param malloc_type : type of memory allacation : MA_MALLOC, MA_PAGELOCK,
   * MA_ZERO, MA_PORT, MA_WC
   */
  init(dims_data, 0L, malloc_type, nb_streams);
}
template CarmaHostObj<float>::CarmaHostObj(const int64_t *dims_data,
                                               MemAlloc malloc_type,
                                               int32_t nb_streams);
template CarmaHostObj<double>::CarmaHostObj(const int64_t *dims_data,
                                                MemAlloc malloc_type,
                                                int32_t nb_streams);
template CarmaHostObj<int32_t>::CarmaHostObj(const int64_t *dims_data,
                                             MemAlloc malloc_type,
                                             int32_t nb_streams);
template CarmaHostObj<uint32_t>::CarmaHostObj(const int64_t *dims_data,
                                                      MemAlloc malloc_type,
                                                      int32_t nb_streams);
template CarmaHostObj<cuFloatComplex>::CarmaHostObj(const int64_t *dims_data,
                                                        MemAlloc malloc_type,
                                                        int32_t nb_streams);
template CarmaHostObj<cuDoubleComplex>::CarmaHostObj(const int64_t *dims_data,
                                                         MemAlloc malloc_type,
                                                         int32_t nb_streams);

template <class T_data>
CarmaHostObj<T_data>::CarmaHostObj(const CarmaHostObj<T_data> *src,
                                       int32_t nb_streams) {
  /** \brief CarmaHostObj creator.
   * \param src : CarmaHostObj to copy
   */
  init(src->dims_data, src->h_data, src->malloc_type, nb_streams);
}
template CarmaHostObj<float>::CarmaHostObj(const CarmaHostObj<float> *src,
                                               int32_t nb_streams);
template CarmaHostObj<double>::CarmaHostObj(
    const CarmaHostObj<double> *src, int32_t nb_streams);
template CarmaHostObj<int32_t>::CarmaHostObj(const CarmaHostObj<int32_t> *src,
                                             int32_t nb_streams);
template CarmaHostObj<uint32_t>::CarmaHostObj(
    const CarmaHostObj<uint32_t> *src, int32_t nb_streams);
template CarmaHostObj<cuFloatComplex>::CarmaHostObj(
    const CarmaHostObj<cuFloatComplex> *src, int32_t nb_streams);
template CarmaHostObj<cuDoubleComplex>::CarmaHostObj(
    const CarmaHostObj<cuDoubleComplex> *src, int32_t nb_streams);

template <class T_data>
CarmaHostObj<T_data>::CarmaHostObj(const CarmaHostObj<T_data> *src,
                                       MemAlloc malloc_type, int32_t nb_streams) {
  /** \brief CarmaHostObj creator.
   * \param src : CarmaHostObj to copy
   * \param malloc_type : type of memory allacation : MA_MALLOC, MA_PAGELOCK,
   * MA_ZERO, MA_PORT, MA_WC
   */
  init(src->dims_data, src->h_data, malloc_type, nb_streams);
}
template CarmaHostObj<float>::CarmaHostObj(const CarmaHostObj<float> *src,
                                               MemAlloc malloc_type,
                                               int32_t nb_streams);
template CarmaHostObj<double>::CarmaHostObj(
    const CarmaHostObj<double> *src, MemAlloc malloc_type, int32_t nb_streams);
template CarmaHostObj<int32_t>::CarmaHostObj(const CarmaHostObj<int32_t> *src,
                                             MemAlloc malloc_type,
                                             int32_t nb_streams);
template CarmaHostObj<uint32_t>::CarmaHostObj(
    const CarmaHostObj<uint32_t> *src, MemAlloc malloc_type,
    int32_t nb_streams);
template CarmaHostObj<cuFloatComplex>::CarmaHostObj(
    const CarmaHostObj<cuFloatComplex> *src, MemAlloc malloc_type,
    int32_t nb_streams);
template CarmaHostObj<cuDoubleComplex>::CarmaHostObj(
    const CarmaHostObj<cuDoubleComplex> *src, MemAlloc malloc_type,
    int32_t nb_streams);

template <class T_data>
CarmaHostObj<T_data>::CarmaHostObj(const int64_t *dims_data,
                                       const T_data *data, int32_t nb_streams) {
  /** \brief CarmaHostObj creator.
   * \param dims_data : the array size Yorick format : [ndims,dims1,dims2,...]
   * \param data : the array
   */
  init(dims_data, data, MA_MALLOC, nb_streams);
}
template CarmaHostObj<float>::CarmaHostObj(const int64_t *dims_data,
                                               const float *data,
                                               int32_t nb_streams);
template CarmaHostObj<double>::CarmaHostObj(const int64_t *dims_data,
                                                const double *data,
                                                int32_t nb_streams);
template CarmaHostObj<int32_t>::CarmaHostObj(const int64_t *dims_data,
                                             const int32_t *data, int32_t nb_streams);
template CarmaHostObj<uint32_t>::CarmaHostObj(const int64_t *dims_data,
                                                      const uint32_t *data,
                                                      int32_t nb_streams);
template CarmaHostObj<cuFloatComplex>::CarmaHostObj(
    const int64_t *dims_data, const cuFloatComplex *data, int32_t nb_streams);
template CarmaHostObj<cuDoubleComplex>::CarmaHostObj(
    const int64_t *dims_data, const cuDoubleComplex *data, int32_t nb_streams);

template <class T_data>
CarmaHostObj<T_data>::CarmaHostObj(const int64_t *dims_data,
                                       const T_data *data, MemAlloc malloc_type,
                                       int32_t nb_streams) {
  /** \brief CarmaHostObj creator.
   * \param dims_data : the array size Yorick format : [ndims,dims1,dims2,...]
   * \param data : the array
   * \param malloc_type : type of memory allacation : MA_MALLOC, MA_PAGELOCK,
   * MA_ZERO, MA_PORT, MA_WC
   */
  init(dims_data, data, malloc_type, nb_streams);
}
template CarmaHostObj<float>::CarmaHostObj(const int64_t *dims_data,
                                               const float *data,
                                               MemAlloc malloc_type,
                                               int32_t nb_streams);
template CarmaHostObj<double>::CarmaHostObj(const int64_t *dims_data,
                                                const double *data,
                                                MemAlloc malloc_type,
                                                int32_t nb_streams);
template CarmaHostObj<int32_t>::CarmaHostObj(const int64_t *dims_data,
                                             const int32_t *data,
                                             MemAlloc malloc_type,
                                             int32_t nb_streams);
template CarmaHostObj<uint32_t>::CarmaHostObj(const int64_t *dims_data,
                                                      const uint32_t *data,
                                                      MemAlloc malloc_type,
                                                      int32_t nb_streams);
template CarmaHostObj<cuFloatComplex>::CarmaHostObj(
    const int64_t *dims_data, const cuFloatComplex *data, MemAlloc malloc_type,
    int32_t nb_streams);
template CarmaHostObj<cuDoubleComplex>::CarmaHostObj(
    const int64_t *dims_data, const cuDoubleComplex *data, MemAlloc malloc_type,
    int32_t nb_streams);

template <class T_data>
void CarmaHostObj<T_data>::init(const int64_t *dims_data, const T_data *data,
                                  MemAlloc malloc_type, int32_t nb_streams) {
  const int64_t size_data = dims_data[0] + 1;
  this->dims_data = new int64_t[size_data];
  std::copy(dims_data, dims_data + size_data, this->dims_data);

  this->nb_elem = dims_data[1];
  for (int32_t i = 2; i < size_data; i++) this->nb_elem *= dims_data[i];

  this->malloc_type = malloc_type;
  if (malloc_type == MA_MALLOC) {
    h_data = new T_data[this->nb_elem];
  } else if (malloc_type == MA_PAGELOCK) {
    carma_safe_call(cudaHostAlloc((void **)&(this->h_data),
                                sizeof(T_data) * this->nb_elem,
                                cudaHostAllocDefault));
  } else if (malloc_type == MA_ZEROCPY) {
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    if (!prop.canMapHostMemory) {
      DEBUG_TRACE("Can't map host memory");
      throw "Can't map host memory\n";
    }
    cudaSetDeviceFlags(cudaDeviceMapHost);
    carma_safe_call(
        cudaHostAlloc((void **)&(this->h_data), sizeof(T_data) * this->nb_elem,
                      cudaHostAllocWriteCombined | cudaHostAllocMapped));
  } else if (malloc_type == MA_PORTABLE) {
    carma_safe_call(
        cudaHostAlloc((void **)&(this->h_data), sizeof(T_data) * this->nb_elem,
                      cudaHostAllocWriteCombined | cudaHostAllocMapped |
                          cudaHostAllocPortable));
  } else if (malloc_type == MA_WRICOMB) {
    carma_safe_call(cudaHostAlloc((void **)&(this->h_data),
                                sizeof(T_data) * this->nb_elem,
                                cudaHostAllocWriteCombined));
  } else if (malloc_type == MA_GENEPIN) {
    cudaSetDeviceFlags(cudaDeviceBlockingSync | cudaDeviceMapHost);
    this->data_UA = (T_data *)mmap(
        NULL, (sizeof(T_data) * this->nb_elem + MEMORY_ALIGNMENT),
        PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANON, -1, 0);
    this->h_data = (T_data *)ALIGN_UP(data_UA, MEMORY_ALIGNMENT);
    carma_safe_call(cudaHostRegister(h_data, sizeof(T_data) * this->nb_elem,
                                   cudaHostRegisterMapped));
  } else
    throw "Error : type of malloc unknown";

  streams = new CarmaStreams();
  if (malloc_type != MA_MALLOC) {
    this->add_stream(nb_streams);
  }

#if DEBUG
  printf("CARMA Host Object created @ 0x%p with %s memory\n", this,
         get_mem_alloc().c_str());
#endif
  if (data == 0L) {
    if (malloc_type == MA_MALLOC) {
      memset(this->h_data, 0, sizeof(T_data) * this->nb_elem);
    } else {
      // carma_safe_call(cudaMemset((T_data *)this->h_data, 0,
      // sizeof(T_data)*this->nb_elem));
      memset(this->h_data, 0, sizeof(T_data) * this->nb_elem);
    }
  } else
    this->fill_from(data);
}

template <class T_data>
void CarmaHostObj<T_data>::get_devpntr(void **pntr_dev) {
  /** \brief retreive device pointer from a host object
   */
  // carma_safe_call(cudaHostGetDevicePointer((void **)&d_a, (void
  // *)(this->d_data), this->malloc_type));
  carma_safe_call(cudaHostGetDevicePointer(pntr_dev, (void *)(this->h_data), 0));
}

/*
 template<class T_data>
 CarmaHostObj<T_data>& CarmaHostObj<T_data>::operator= (const
 CarmaHostObj<T_data>& obj){

 if(this->d_data!=0L) carma_safe_call( cudaFree(this->d_data) );
 if(this->dims_data!=0L) delete(this->dims_data);
 CarmaHostObj<T_data> *new_obj = new CarmaHostObj(obj);
 return new_obj;
 }
 */

template <class T_data>
CarmaHostObj<T_data>::~CarmaHostObj() {
  /** \brief CarmaHostObj destructor.
   */

  if (this->h_data != 0L) {
    if (malloc_type == MA_MALLOC)
      delete[] this->h_data;
    else if (malloc_type == MA_GENEPIN) {
      carma_safe_call(cudaHostUnregister(this->h_data));
      munmap(this->data_UA, sizeof(T_data) * this->nb_elem);
    } else
      carma_safe_call(cudaFreeHost(this->h_data));
    this->h_data = 0L;
  }

  delete this->streams;

#if DEBUG
  printf("CARMA Host Object deleted @ 0x%p\n", this);
#endif
}
template CarmaHostObj<float>::~CarmaHostObj();
template CarmaHostObj<double>::~CarmaHostObj();
template CarmaHostObj<int32_t>::~CarmaHostObj();
template CarmaHostObj<uint32_t>::~CarmaHostObj();
template CarmaHostObj<cuFloatComplex>::~CarmaHostObj();
template CarmaHostObj<cuDoubleComplex>::~CarmaHostObj();

template <class T_data>
int32_t CarmaHostObj<T_data>::get_nb_streams() {
  /** \brief get the number of streams attached to the host object
   */
  return streams->get_nb_streams();
}

template int32_t CarmaHostObj<float>::get_nb_streams();
template int32_t CarmaHostObj<double>::get_nb_streams();
template int32_t CarmaHostObj<int32_t>::get_nb_streams();
template int32_t CarmaHostObj<uint32_t>::get_nb_streams();
template int32_t CarmaHostObj<cuFloatComplex>::get_nb_streams();
template int32_t CarmaHostObj<cuDoubleComplex>::get_nb_streams();

template <class T_data>
int32_t CarmaHostObj<T_data>::add_stream() {
  if (this->malloc_type != MA_MALLOC)
    this->streams->add_stream();
  else
    std::cerr << "malloc Carma Host Object can not be streamed" << std::endl;
  return this->streams->get_nb_streams();
}

template <class T_data>
int32_t CarmaHostObj<T_data>::add_stream(int32_t nb) {
  if (this->malloc_type != MA_MALLOC)
    this->streams->add_stream(nb);
  else
    std::cerr << "malloc Carma Host Object can not be streamed" << std::endl;
  return this->streams->get_nb_streams();
}

template <class T_data>
cudaStream_t CarmaHostObj<T_data>::get_cuda_stream(int32_t stream) {
  return this->streams->get_stream(stream);
}

template cudaStream_t CarmaHostObj<float>::get_cuda_stream(int32_t stream);
template cudaStream_t CarmaHostObj<double>::get_cuda_stream(int32_t stream);
template cudaStream_t CarmaHostObj<int32_t>::get_cuda_stream(int32_t stream);
template cudaStream_t CarmaHostObj<uint32_t>::get_cuda_stream(
    int32_t stream);
template cudaStream_t CarmaHostObj<cuFloatComplex>::get_cuda_stream(
    int32_t stream);
template cudaStream_t CarmaHostObj<cuDoubleComplex>::get_cuda_stream(
    int32_t stream);

template <class T_data>
int32_t CarmaHostObj<T_data>::del_stream() {
  this->streams->del_stream();
  return this->streams->get_nb_streams();
}

template <class T_data>
int32_t CarmaHostObj<T_data>::del_stream(int32_t nb) {
  this->streams->del_stream(nb);
  return this->streams->get_nb_streams();
}

template <class T_data>
int32_t CarmaHostObj<T_data>::wait_stream(int32_t stream) {
  this->streams->wait_stream(stream);
  return EXIT_SUCCESS;
}

template <class T_data>
int32_t CarmaHostObj<T_data>::wait_all_streams() {
  this->streams->wait_all_streams();
  return EXIT_SUCCESS;
}
template int32_t CarmaHostObj<float>::wait_all_streams();

/*
 template<class T_data>
 T_data* CarmaHostObj<T_data>::get_data(){
 * \brief get_data data transfer.
 * \return data : pointer on the data
 *

 return h_data;
 }
 */
template <class T_data>
int32_t CarmaHostObj<T_data>::fill(T_data value) {
  std::fill(this->h_data, this->h_data + this->get_nb_elements(), value);

  return EXIT_SUCCESS;
}

template int32_t CarmaHostObj<float>::fill(float value);
template int32_t CarmaHostObj<double>::fill(double value);
template int32_t CarmaHostObj<int32_t>::fill(int32_t value);
template int32_t CarmaHostObj<uint32_t>::fill(uint32_t value);
template int32_t CarmaHostObj<cuFloatComplex>::fill(cuFloatComplex value);
template int32_t CarmaHostObj<cuDoubleComplex>::fill(cuDoubleComplex value);

template <class T_data>
int32_t CarmaHostObj<T_data>::fill_from(const T_data *data) {
  if (get_nb_streams() > 1) {
    int32_t nstreams = get_nb_streams();
    for (int32_t i = 0; i < nstreams; i++) {
      carma_safe_call(cudaMemcpyAsync(&(h_data[i * nb_elem / nstreams]),
                                    &(data[i * nb_elem / nstreams]),
                                    sizeof(T_data) * nb_elem / nstreams,
                                    cudaMemcpyHostToHost, get_cuda_stream(i)));
    }
  } else if (get_nb_streams() == 1) {
    carma_safe_call(cudaMemcpyAsync(h_data, data, nb_elem * sizeof(T_data),
                                  cudaMemcpyHostToHost, get_cuda_stream(0)));
  } else {
    if (malloc_type == MA_MALLOC) {
      std::copy(data, data + nb_elem, h_data);
    } else {
      carma_safe_call(cudaMemcpy(h_data, data, nb_elem * sizeof(T_data),
                               cudaMemcpyHostToHost));
    }
  }
  return EXIT_SUCCESS;
}
template int32_t CarmaHostObj<float>::fill_from(const float *data);
template int32_t CarmaHostObj<double>::fill_from(const double *data);
template int32_t CarmaHostObj<int32_t>::fill_from(const int32_t *data);
template int32_t CarmaHostObj<uint32_t>::fill_from(const uint32_t *data);
template int32_t CarmaHostObj<cuFloatComplex>::fill_from(
    const cuFloatComplex *data);
template int32_t CarmaHostObj<cuDoubleComplex>::fill_from(
    const cuDoubleComplex *data);

template <class T_data>
int32_t CarmaHostObj<T_data>::fill_into(T_data *data) {
  if (get_nb_streams() > 1) {
    int32_t nstreams = get_nb_streams();
    for (int32_t i = 0; i < nstreams; i++) {
      carma_safe_call(cudaMemcpyAsync(&(data[i * nb_elem / nstreams]),
                                    &(h_data[i * nb_elem / nstreams]),
                                    sizeof(T_data) * nb_elem / nstreams,
                                    cudaMemcpyHostToHost, get_cuda_stream(i)));
    }
  } else if (get_nb_streams() == 1) {
    carma_safe_call(cudaMemcpyAsync(data, h_data, nb_elem * sizeof(T_data),
                                  cudaMemcpyHostToHost, get_cuda_stream(0)));
  } else {
    if (malloc_type == MA_MALLOC) {
      std::copy(h_data, h_data + nb_elem, data);
    } else {
      carma_safe_call(cudaMemcpy(data, h_data, nb_elem * sizeof(T_data),
                               cudaMemcpyHostToHost));
    }
  }
  return EXIT_SUCCESS;
}
template int32_t CarmaHostObj<float>::fill_into(float *data);
template int32_t CarmaHostObj<double>::fill_into(double *data);
template int32_t CarmaHostObj<int32_t>::fill_into(int32_t *data);
template int32_t CarmaHostObj<uint32_t>::fill_into(uint32_t *data);
template int32_t CarmaHostObj<cuFloatComplex>::fill_into(cuFloatComplex *data);
template int32_t CarmaHostObj<cuDoubleComplex>::fill_into(cuDoubleComplex *data);

template <class T_data>
int32_t CarmaHostObj<T_data>::cpy_obj(CarmaObj<T_data> *carma_obj,
                                    cudaMemcpyKind flag) {
  if (nb_elem != carma_obj->get_nb_elements()) {
    std::cerr << "***** ERROR (" << __FILE__ << "@" << __LINE__
              << ") : two objects do not have the same size *****" << std::endl;
    return EXIT_FAILURE;
  }

  T_data *data_src, *data_dst;
  if (flag == cudaMemcpyHostToDevice) {
    data_dst = *carma_obj;
    data_src = h_data;
  } else if (flag == cudaMemcpyDeviceToHost) {
    data_dst = h_data;
    data_src = *carma_obj;
  } else {
    std::cerr << "***** ERROR (" << __FILE__ << "@" << __LINE__
              << ") : wrong flag *****" << std::endl;
    return EXIT_FAILURE;
  }

  if (get_nb_streams() > 1) {
    int32_t nstreams = get_nb_streams();
    for (int32_t i = 0; i < nstreams; i++) {
      carma_safe_call(cudaMemcpyAsync(&(data_dst[i * nb_elem / nstreams]),
                                    &(data_src[i * nb_elem / nstreams]),
                                    sizeof(T_data) * nb_elem / nstreams, flag,
                                    get_cuda_stream(i)));
    }
  } else if (get_nb_streams() == 1) {
    carma_safe_call(cudaMemcpyAsync(data_dst, data_src, sizeof(T_data) * nb_elem,
                                  flag, get_cuda_stream(0)));
  } else {
    carma_safe_call(
        cudaMemcpy(data_dst, data_src, nb_elem * sizeof(T_data), flag));
  }
  return EXIT_SUCCESS;
}
template int32_t CarmaHostObj<float>::cpy_obj(CarmaObj<float> *data,
                                            cudaMemcpyKind flag);
template int32_t CarmaHostObj<double>::cpy_obj(CarmaObj<double> *data,
                                             cudaMemcpyKind flag);
template int32_t CarmaHostObj<int32_t>::cpy_obj(CarmaObj<int32_t> *data,
                                          cudaMemcpyKind flag);
template int32_t CarmaHostObj<uint32_t>::cpy_obj(
    CarmaObj<uint32_t> *data, cudaMemcpyKind flag);
template int32_t CarmaHostObj<cuFloatComplex>::cpy_obj(
    CarmaObj<cuFloatComplex> *data, cudaMemcpyKind flag);
template int32_t CarmaHostObj<cuDoubleComplex>::cpy_obj(
    CarmaObj<cuDoubleComplex> *data, cudaMemcpyKind flag);

template <class T_data>
int32_t CarmaHostObj<T_data>::cpy_obj(CarmaObj<T_data> *carma_obj,
                                    cudaMemcpyKind flag, uint32_t stream) {
  uint32_t nb_streams = this->streams->get_nb_streams();
  if (stream >= nb_streams) return EXIT_FAILURE;

  CarmaStreams streams_tmp = *(this->streams);
  if (flag == cudaMemcpyHostToDevice) {
    carma_safe_call(cudaMemcpyAsync(*carma_obj, h_data,
                                  carma_obj->get_nb_elements() * sizeof(T_data), flag,
                                  get_cuda_stream(stream)));
  } else if (flag == cudaMemcpyDeviceToHost) {
    carma_safe_call(cudaMemcpyAsync(h_data, *carma_obj,
                                  carma_obj->get_nb_elements() * sizeof(T_data), flag,
                                  get_cuda_stream(stream)));
  } else
    return EXIT_FAILURE;
  return EXIT_SUCCESS;
}
template int32_t CarmaHostObj<float>::cpy_obj(CarmaObj<float> *data,
                                            cudaMemcpyKind flag,
                                            uint32_t stream);
template int32_t CarmaHostObj<double>::cpy_obj(CarmaObj<double> *data,
                                             cudaMemcpyKind flag,
                                             uint32_t stream);
template int32_t CarmaHostObj<int32_t>::cpy_obj(CarmaObj<int32_t> *data,
                                          cudaMemcpyKind flag,
                                          uint32_t stream);
template int32_t CarmaHostObj<uint32_t>::cpy_obj(
    CarmaObj<uint32_t> *data, cudaMemcpyKind flag, uint32_t stream);
template int32_t CarmaHostObj<cuFloatComplex>::cpy_obj(
    CarmaObj<cuFloatComplex> *data, cudaMemcpyKind flag, uint32_t stream);
template int32_t CarmaHostObj<cuDoubleComplex>::cpy_obj(
    CarmaObj<cuDoubleComplex> *data, cudaMemcpyKind flag, uint32_t stream);
