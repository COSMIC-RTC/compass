// This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
//
// COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
// General Public License as published by the Free Software Foundation, either version 3 of the 
// License, or any later version.
//
// COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
// See the GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License along with COMPASS. 
// If not, see <https://www.gnu.org/licenses/>
//
//  Copyright (C) 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>

//! \file      carma_obj.hpp
//! \ingroup   libcarma
//! \class     CarmaObj
//! \brief     this class provides wrappers to the generic carma object
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#ifndef _CARMA_OBJ_H_
#define _CARMA_OBJ_H_

#include <carma_context.hpp>
#include <carma_cublas.hpp>
#include <carma_streams.hpp>
#include <carma_utils.hpp>
#include <curand.h>
#include <curand_kernel.h>
#include <iostream>
#include <type_traits>
#include <typeinfo>  // operator typeid

/*
 create a memory object
 void *memory
 int32_t  nb of reference

 create a class which contains :
 - d_data
 - ndims
 - dims
 - strides
 - type

 new()

 new(existing)


 and then
 modify CarmaObj so that it is :
 an object of the previous class
 all the methods of a CarmaObj

 */

#define BLOCK_SZ 16

enum MemType {
  MT_DEVICE,
  MT_DARRAY,
  MT_HOST,
  MT_PAGELOCK,
  MT_ZEROCPY,
  MT_PORTABLE,
  MT_WRICOMB,
  MT_GENEPIN
};
// should add texture ?

template <class T_data>
class CarmaData {
 protected:
  T_data *d_data;       ///< Pointer to data
  int32_t ndims;            ///< Number of dimensions
  int32_t nb_elem;          ///< Number of elements
  int64_t *dims_data;      ///< Dimensions
  int32_t *strides;         ///< Strides for each dimension
  MemType malloc_type;  ///< type of alloc

 public:
  T_data *get_data() { return d_data; }
  int32_t get_ndims() { return ndims; }
  int32_t get_nb_elem() { return nb_elem; }
  const int64_t *get_dims_data() { return dims_data; }
  int64_t get_dims_data(int32_t i) { return dims_data[i]; }
  int32_t *get_strides() { return strides; }
  int32_t get_strides(int32_t i) { return strides[i]; }
  MemType get_malloc_type() { return malloc_type; }
};

template <class T_data>
class CarmaHostObj;

template <class T_data>
class CarmaObj {
 protected:
  T_data *d_data = nullptr;  ///< Input data  => change to vector
  std::vector<T_data> h_data;
  T_data *o_data = nullptr;        ///< optional data (used for scan / reduction)
  T_data *cub_data = nullptr;      ///< optional data (used for scan / reduction)
  size_t cub_data_size = 0;  // optionnal for reduction
  int32_t ndim = 0;
  int64_t *dims_data = nullptr;  ///< dimensions of the array
  int32_t nb_elem = 0;      ///< number of elements in the array
  int32_t device = -1;       ///< device where the CarmaObj is allocate
  CarmaContext *current_context;

  curandGenerator_t gen;
  curandState *d_states;

  int32_t nb_threads;
  int32_t nb_blocks;

  bool keys_only;      //< optional flag (used for sort)
  bool owner = true;  // Flag if d_data is created inside the CarmaObj

  uint32_t *values;  ///< optional data (used for sort)
  size_t *d_num_valid;    ///< used for compact

  cufftHandle plan;  ///< FFT plan
  cufftType type_plan;   ///< FFT plan type

  CarmaStreams *streams;

  void init(CarmaContext *current_context, const int64_t *dims_data,
            const T_data *data, bool fromHost, int32_t nb_streams);

 public:
  CarmaObj();
  CarmaObj(const CarmaObj<T_data> *obj);
  CarmaObj(CarmaContext *current_context, const int64_t *dims_data);
  CarmaObj(CarmaContext *current_context, const std::vector<int64_t> &dims);
  CarmaObj(CarmaContext *current_context, const CarmaObj<T_data> *obj);
  CarmaObj(CarmaContext *current_context, const int64_t *dims_data,
            const T_data *data);
  CarmaObj(CarmaContext *current_context, const int64_t *dims_data,
            int32_t nb_streams);
  CarmaObj(CarmaContext *current_context, const CarmaObj<T_data> *obj,
            int32_t nb_streams);
  CarmaObj(CarmaContext *current_context, const int64_t *dims_data,
            const T_data *data, int32_t nb_streams);
  CarmaObj(const CarmaObj &)=delete;
  ~CarmaObj();

  void sync_h_data() {
    if (h_data.empty()) h_data = std::vector<T_data>(nb_elem);
    device2host(h_data.data());
  }

  T_data *get_h_data() { return h_data.data(); }

  int32_t get_nb_streams() const {
    /** \brief get the number of streams attached to the host object
     */
    return streams->get_nb_streams();
  }
  int32_t add_stream() {
    this->streams->add_stream();
    return this->streams->get_nb_streams();
  }
  int32_t add_stream(int32_t nb) {
    this->streams->add_stream(nb);
    return this->streams->get_nb_streams();
  }
  int32_t del_stream() {
    this->streams->del_stream();
    return this->streams->get_nb_streams();
  }
  int32_t del_stream(int32_t nb) {
    this->streams->del_stream(nb);
    return this->streams->get_nb_streams();
  }
  cudaStream_t get_cuda_stream(int32_t stream) {
    return this->streams->get_stream(stream);
  }
  int32_t wait_stream(int32_t stream) {
    this->streams->wait_stream(stream);
    return EXIT_SUCCESS;
  }
  int32_t wait_all_streams() {
    this->streams->wait_all_streams();
    return EXIT_SUCCESS;
  }
  void swap_ptr(T_data *ptr) {
    dealloc();
    d_data = ptr;
    owner = false;
  }

  void dealloc() {
    if (owner && d_data) cudaFree(d_data);
  }

  /**< General Utilities */
  operator T_data *() { return d_data; }

  std::string to_string() {
    std::ostringstream stream;
    stream << *this;
    return stream.str();
  }

  operator std::string() { return this->to_string(); }
  // inline char const *c_str() { return this->to_string().c_str(); }
  const T_data operator[](int32_t index) const {
    T_data tmp_float;
    carma_safe_call(cudaMemcpy(&tmp_float, &d_data[index], sizeof(T_data),
                             cudaMemcpyDeviceToHost));
    return tmp_float;
  }
  T_data *get_data() { return d_data; }
  T_data *get_data_at(int32_t index) { return &d_data[index]; }
  T_data *get_o_data() { return o_data; }
  const T_data get_o_data_value() const {
    T_data tmp_float;
    carma_safe_call(
        cudaMemcpy(&tmp_float, o_data, sizeof(T_data), cudaMemcpyDeviceToHost));
    return tmp_float;
  }
  const int64_t *get_dims() { return dims_data; }
  int64_t get_dims(int32_t i) { return dims_data[i]; }
  int32_t get_nb_elements() { return nb_elem; }
  CarmaContext *get_context() { return current_context; }

  int32_t get_device() { return device; }

  bool is_rng_init() { return (gen != NULL); }

  /**< Memory transfers both ways */
  template <typename T_dest>
  int32_t host2device(const T_dest *data);
  template <typename T_dest>
  int32_t device2host(T_dest *data);

  int32_t host2device_async(const T_data *data, cudaStream_t stream);
  int32_t device2host_async(T_data *data, cudaStream_t stream);
  int32_t device2host_opt(T_data *data);
  int32_t host2device_vect(const T_data *data, int32_t incx, int32_t incy);
  int32_t device2host_vect(T_data *data, int32_t incx, int32_t incy);
  int32_t host2device_mat(const T_data *data, int32_t lda, int32_t ldb);
  int32_t device2host_mat(T_data *data, int32_t lda, int32_t ldb);

  int32_t copy_into(T_data *data, int32_t nb_elem);
  int32_t copy_from(const T_data *data, int32_t nb_elem);
  int32_t copy_from_async(const T_data *data, int32_t nb_elem, cudaStream_t stream);

  inline int32_t reset() {
    return cudaMemset(this->d_data, 0, this->nb_elem * sizeof(T_data));
  }

  inline int32_t reset(cudaStream_t stream) {
    return cudaMemsetAsync(this->d_data, 0, this->nb_elem * sizeof(T_data), stream);
  }

  inline int32_t memset(T_data value) {
    return fill_array_with_value(
        this->d_data, value, this->nb_elem,
        this->current_context->get_device(this->device));
  }
  cufftHandle *get_plan() { return &plan; }
  ///< FFT plan
  cufftType get_type_plan() { return type_plan; }
  ///< FFT plan type

  uint32_t *get_values() { return values; }
  ///< optional data (used for sort)

  /**< sum */
  T_data sum();
  void init_reduceCub();
  void reduceCub(cudaStream_t stream);
  void reduceCub() {reduceCub(0);};

  void clip(T_data min, T_data max, cudaStream_t stream);
  void clip(T_data min, T_data max) {clip(min, max, 0);};

  /**< transpose */
  int32_t transpose(CarmaObj<T_data> *source);
  // CarmaObj<T_data>& operator= (const CarmaObj<T_data>& obj);

  /*
   *  ____  _        _    ____  _
   * | __ )| |      / \  / ___|/ |
   * |  _ \| |     / _ \ \___ \| |
   * | |_) | |___ / ___ \ ___) | |
   * |____/|_____/_/   \_\____/|_|
   *
   */

  int32_t aimax(int32_t incx);
  int32_t aimin(int32_t incx);
  T_data asum(int32_t incx);
  T_data nrm2(int32_t incx);
  T_data dot(CarmaObj<T_data> *source, int32_t incx, int32_t incy);
  void scale(T_data alpha, int32_t incx);
  void swap(CarmaObj<T_data> *source, int32_t incx, int32_t incy);
  void copy(CarmaObj<T_data> *source, int32_t incx, int32_t incy);
  void axpy(T_data alpha, CarmaObj<T_data> *source, int32_t incx, int32_t incy,
            int32_t offset = 0);
  void rot(CarmaObj<T_data> *source, int32_t incx, int32_t incy, T_data sc, T_data ss);

  /*
   *  ____  _        _    ____ ____
   * | __ )| |      / \  / ___|___ \
   * |  _ \| |     / _ \ \___ \ __) |
   * | |_) | |___ / ___ \ ___) / __/
   * |____/|_____/_/   \_\____/_____|
   *
   */

  void gemv(char trans, T_data alpha, CarmaObj<T_data> *matA, int32_t lda,
            CarmaObj<T_data> *vectx, int32_t incx, T_data beta, int32_t incy);
  void ger(T_data alpha, CarmaObj<T_data> *vectx, int32_t incx,
           CarmaObj<T_data> *vecty, int32_t incy, int32_t lda);
  void symv(char uplo, T_data alpha, CarmaObj<T_data> *matA, int32_t lda,
            CarmaObj<T_data> *vectx, int32_t incx, T_data beta, int32_t incy);

  /*
   *  ____  _        _    ____ _____
   * | __ )| |      / \  / ___|___ /
   * |  _ \| |     / _ \ \___ \ |_ \
   * | |_) | |___ / ___ \ ___) |__) |
   * |____/|_____/_/   \_\____/____/
   *
   */

  void gemm(char transa, char transb, T_data alpha, CarmaObj<T_data> *matA,
            int32_t lda, CarmaObj<T_data> *matB, int32_t ldb, T_data beta, int32_t ldc);
  void symm(char side, char uplo, T_data alpha, CarmaObj<T_data> *matA,
            int32_t lda, CarmaObj<T_data> *matB, int32_t ldb, T_data beta, int32_t ldc);
  void syrk(char uplo, char transa, T_data alpha, CarmaObj<T_data> *matA,
            int32_t lda, T_data beta, int32_t ldc);
  void syrkx(char uplo, char transa, T_data alpha, CarmaObj<T_data> *matA,
             int32_t lda, CarmaObj<T_data> *matB, int32_t ldb, T_data beta, int32_t ldc);
  void geam(char transa, char transb, T_data alpha, CarmaObj<T_data> *matA,
            int32_t lda, T_data beta, CarmaObj<T_data> *matB, int32_t ldb, int32_t ldc);
  void dgmm(char side, CarmaObj<T_data> *matA, int32_t lda,
            CarmaObj<T_data> *vectx, int32_t incx, int32_t ldc);

  /**< Curand */
  int32_t init_prng();
  int32_t init_prng(int64_t seed);
  int32_t destroy_prng();
  int32_t prng(T_data *output, char gtype, float alpha, float beta);
  int32_t prng(T_data *output, char gtype, float alpha);
  int32_t prng(char gtype, float alpha, float beta);
  int32_t prng(char gtype, float alpha);
  int32_t prng(char gtype);

  int32_t prng_montagn(float init_montagn);

  int32_t init_prng_host(int32_t seed);
  int32_t prng_host(char gtype);
  int32_t prng_host(char gtype, T_data stddev);
  int32_t prng_host(char gtype, T_data stddev, T_data alpha);
  int32_t destroy_prng_host();
};

typedef CarmaObj<int32_t> CarmaObjI;
typedef CarmaObj<uint32_t> CarmaObjUI;
typedef CarmaObj<uint16_t> CarmaObjUSI;
typedef CarmaObj<float> CarmaObjS;
typedef CarmaObj<double> CarmaObjD;
typedef CarmaObj<float2> CarmaObjS2;
typedef CarmaObj<double2> CarmaObjD2;
typedef CarmaObj<cuFloatComplex> CarmaObjC;
typedef CarmaObj<cuDoubleComplex> CarmaObjZ;
// typedef CarmaObj<tuple_t<float>> CarmaObjTF;

template <class T_data>
std::ostream &operator<<(std::ostream &os, CarmaObj<T_data> &obj) {
  os << "-----------------------" << std::endl;
  os << "CarmaObj<" << typeid(T_data).name() << "> object on GPU"
     << obj.get_device() << std::endl;
  int64_t ndims = obj.get_dims(0);
  os << "ndims = " << ndims << std::endl;
  for (int64_t dim = 0; dim < ndims; dim++) {
    os << "dim[" << dim << "] = " << obj.get_dims(dim + 1) << std::endl;
  }
  os << "nbElem = " << obj.get_nb_elements() << std::endl;
  os << "sizeof(" << typeid(T_data).name() << ") = " << sizeof(T_data)
     << std::endl;
  os << "-----------------------" << std::endl;
  return os;
}

// CU functions clip
template <class T_data>
void clip_array(T_data *d_data, T_data min, T_data max, int32_t N,
                CarmaDevice *device, cudaStream_t stream);

// CU functions sum
template <class T_data>
void reduce(int32_t size, int32_t threads, int32_t blocks, T_data *d_idata,
            T_data *d_odata);
template <class T_data>
T_data reduce(T_data *data, int32_t N);

template <class T_data>
void init_reduceCubCU(T_data *&cub_data, size_t &cub_data_size, T_data *data,
                      T_data *&o_data, int32_t N);
template <class T_data>
void reduceCubCU(T_data *cub_data, size_t cub_data_size, T_data *data,
                 T_data *o_data, int32_t N, cudaStream_t stream=0);

// CU functions transpose
template <class T_data>
int32_t transposeCU(T_data *d_idata, T_data *d_odata, int64_t N1, int64_t N2);

// CU functions generic
template <class T_data>
int32_t launch_generic1d(T_data *d_idata, T_data *d_odata, int32_t N,
                     CarmaDevice *device);
template <class T_data>
int32_t launch_generic2d(T_data *d_odata, T_data *d_idata, int32_t N1, int32_t N2);

// CU functions curand
int32_t carma_prng_init(int32_t *seed, const int32_t nb_threads, const int32_t nb_blocks,
                    curandState *state);
template <class T>
int32_t carma_prng_cu(T *results, const int32_t nb_threads, const int32_t nb_blocks,
                  curandState *state, char gtype, int32_t n, float alpha,
                  float beta);
template <class T>
int32_t carma_curand_montagn(curandState *state, T *d_odata, int32_t N,
                         CarmaDevice *device);

// CU functions fft
template <class T_in, class T_out>
cufftType carma_select_plan();
template <class T_in, class T_out>
void carma_initfft(const int64_t *dims_data, cufftHandle *plan, cufftType type_plan);
template <class T_in, class T_out>
int32_t CarmaFFT(T_in *input, T_out *output, int32_t dir, cufftHandle plan);

// CU functions generic
template <class T_data>
int32_t fillindex(T_data *d_odata, T_data *d_idata, int32_t *indx, int32_t N,
              CarmaDevice *device);
template <class T_data>
int32_t fillvalues(T_data *d_odata, T_data *val, int32_t N, CarmaDevice *device, cudaStream_t stream=0);
template <class T>
int32_t getarray2d(T *d_odata, T *d_idata, int32_t x0, int32_t Ncol, int32_t NC, int32_t N,
               CarmaDevice *device);
template <class T>
int32_t fillarray2d(T *d_odata, T *d_idata, int32_t x0, int32_t Ncol, int32_t NC, int32_t N,
                CarmaDevice *device);
template <class T>
int32_t fillarray2d2(T *d_odata, T *d_idata, int32_t x0, int32_t Ncol, int32_t NC, int32_t N,
                 CarmaDevice *device);
template <class T>
int32_t fill_sym_matrix(char src_uplo, T *d_data, int32_t Ncol, int32_t N,
                    CarmaDevice *device);
template <class T>
int32_t carma_plus(T *d_odata, T elpha, int32_t N, CarmaDevice *device);
template <class T>
int32_t carma_plusai(T *d_odata, T *i_data, int32_t i, int32_t sgn, int32_t N,
                 CarmaDevice *device);

// CU functions fftconv
// int32_t fftconv_unpad(float *d_odata, float *d_idata, int32_t fftW, int32_t dataH,
//                   int32_t dataW, int32_t N, int32_t n, int32_t nim);
// int32_t carma_initfftconv(CarmaObjS *data_in, CarmaObjS *kernel_in, CarmaObjS *padded_data,
//                       CarmaObjC *padded_spectrum, int32_t kernelY, int32_t kernelX);
// // CPP functions fftconv
// int32_t carma_fftconv(CarmaObjS *data_out, CarmaObjS *padded_data,
//                   CarmaObjC *padded_spectrum, int32_t kernelY, int32_t kernelX);


/**
 * @brief Kernel to extract a part of the image centred on center_pos
 *
 * @tparam T type of the image items
 * @param d_smallimg extracted small image of size extract_size*extract_size
 * @param d_fullimg full image of size fullimg_size*fullimg_size
 * @param fullimg_size size of the d_fullimg leading dimension
 * @param center_pos position of the center of d_smallimg in d_fullimg
 * @param extract_size size of the d_smallimg leading dimension
 * @param roll get pixels as if d_fullimg need to be roll
 */
template <class T>
int32_t extract(T *d_smallimg, const T *d_fullimg, int32_t fullimg_size, int32_t center_pos,
            int32_t extract_size, bool roll);

#endif  // _CARMA_OBJ_H_
