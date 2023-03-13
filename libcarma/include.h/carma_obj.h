// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      carma_obj.h
//! \ingroup   libcarma
//! \class     CarmaObj
//! \brief     this class provides wrappers to the generic carma object
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.0
//! \date      2022/01/24

#ifndef _CARMA_OBJ_H_
#define _CARMA_OBJ_H_

#include <carma_context.h>
#include <carma_streams.h>
#include <carma_utils.h>
#include <curand.h>
#include <curand_kernel.h>
#include <iostream>
#include <type_traits>
#include <typeinfo>  // operator typeid

/*
 create a memory object
 void *memory
 int  nb of reference

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
  int ndims;            ///< Number of dimensions
  int nb_elem;          ///< Number of elements
  long *dims_data;      ///< Dimensions
  int *strides;         ///< Strides for each dimension
  MemType malloc_type;  ///< type of alloc

 public:
  T_data *get_data() { return d_data; }
  int get_ndims() { return ndims; }
  int get_nb_elem() { return nb_elem; }
  const long *get_dims_data() { return dims_data; }
  long get_dims_data(int i) { return dims_data[i]; }
  int *get_strides() { return strides; }
  int get_strides(int i) { return strides[i]; }
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
  int ndim = 0;
  long *dims_data = nullptr;  ///< dimensions of the array
  int nb_elem = 0;      ///< number of elements in the array
  int device = -1;       ///< device where the CarmaObj is allocate
  CarmaContext *current_context;

  curandGenerator_t gen;
  curandState *d_states;

  int nb_threads;
  int nb_blocks;

  bool keys_only;      //< optional flag (used for sort)
  bool owner = true;  // Flag if d_data is created inside the CarmaObj

  unsigned int *values;  ///< optional data (used for sort)
  size_t *d_num_valid;    ///< used for compact

  cufftHandle plan;  ///< FFT plan
  cufftType type_plan;   ///< FFT plan type

  CarmaStreams *streams;

  void init(CarmaContext *current_context, const long *dims_data,
            const T_data *data, bool fromHost, int nb_streams);

 public:
  CarmaObj();
  CarmaObj(const CarmaObj<T_data> *obj);
  CarmaObj(CarmaContext *current_context, const long *dims_data);
  CarmaObj(CarmaContext *current_context, const std::vector<long> &dims);
  CarmaObj(CarmaContext *current_context, const CarmaObj<T_data> *obj);
  CarmaObj(CarmaContext *current_context, const long *dims_data,
            const T_data *data);
  CarmaObj(CarmaContext *current_context, const long *dims_data,
            int nb_streams);
  CarmaObj(CarmaContext *current_context, const CarmaObj<T_data> *obj,
            int nb_streams);
  CarmaObj(CarmaContext *current_context, const long *dims_data,
            const T_data *data, int nb_streams);
  CarmaObj(const CarmaObj &)=delete;
  ~CarmaObj();

  void sync_h_data() {
    if (h_data.empty()) h_data = std::vector<T_data>(nb_elem);
    device2host(h_data.data());
  }

  T_data *get_h_data() { return h_data.data(); }

  int get_nb_streams() const {
    /** \brief get the number of streams attached to the host object
     */
    return streams->get_nb_streams();
  }
  int add_stream() {
    this->streams->add_stream();
    return this->streams->get_nb_streams();
  }
  int add_stream(int nb) {
    this->streams->add_stream(nb);
    return this->streams->get_nb_streams();
  }
  int del_stream() {
    this->streams->del_stream();
    return this->streams->get_nb_streams();
  }
  int del_stream(int nb) {
    this->streams->del_stream(nb);
    return this->streams->get_nb_streams();
  }
  cudaStream_t get_cuda_stream(int stream) {
    return this->streams->get_stream(stream);
  }
  int wait_stream(int stream) {
    this->streams->wait_stream(stream);
    return EXIT_SUCCESS;
  }
  int wait_all_streams() {
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
  const T_data operator[](int index) const {
    T_data tmp_float;
    carma_safe_call(cudaMemcpy(&tmp_float, &d_data[index], sizeof(T_data),
                             cudaMemcpyDeviceToHost));
    return tmp_float;
  }
  T_data *get_data() { return d_data; }
  T_data *get_data_at(int index) { return &d_data[index]; }
  T_data *get_o_data() { return o_data; }
  const T_data get_o_data_value() const {
    T_data tmp_float;
    carma_safe_call(
        cudaMemcpy(&tmp_float, o_data, sizeof(T_data), cudaMemcpyDeviceToHost));
    return tmp_float;
  }
  const long *get_dims() { return dims_data; }
  long get_dims(int i) { return dims_data[i]; }
  int get_nb_elements() { return nb_elem; }
  CarmaContext *get_context() { return current_context; }

  int get_device() { return device; }

  bool is_rng_init() { return (gen != NULL); }

  /**< Memory transfers both ways */
  template <typename T_dest>
  int host2device(const T_dest *data);
  template <typename T_dest>
  int device2host(T_dest *data);

  int host2device_async(const T_data *data, cudaStream_t stream);
  int device2host_async(T_data *data, cudaStream_t stream);
  int device2host_opt(T_data *data);
  int host2device_vect(const T_data *data, int incx, int incy);
  int device2host_vect(T_data *data, int incx, int incy);
  int host2device_mat(const T_data *data, int lda, int ldb);
  int device2host_mat(T_data *data, int lda, int ldb);

  int copy_into(T_data *data, int nb_elem);
  int copy_from(const T_data *data, int nb_elem);
  int copy_from_async(const T_data *data, int nb_elem, cudaStream_t stream);

#ifdef USE_OCTOPUS
  int copy_into(ipc::Cacao<T_data> *cacaoInterface);
  int copy_from(ipc::Cacao<T_data> *cacaoInterface);
#endif

  inline int reset() {
    return cudaMemset(this->d_data, 0, this->nb_elem * sizeof(T_data));
  }

  inline int reset(cudaStream_t stream) {
    return cudaMemsetAsync(this->d_data, 0, this->nb_elem * sizeof(T_data), stream);
  }

  inline int memset(T_data value) {
    return fill_array_with_value(
        this->d_data, value, this->nb_elem,
        this->current_context->get_device(this->device));
  }
  cufftHandle *get_plan() { return &plan; }
  ///< FFT plan
  cufftType get_type_plan() { return type_plan; }
  ///< FFT plan type

  unsigned int *get_values() { return values; }
  ///< optional data (used for sort)

  /**< sum */
  T_data sum();
  void init_reduceCub();
  void reduceCub(cudaStream_t stream);
  void reduceCub() {reduceCub(0);};

  void clip(T_data min, T_data max, cudaStream_t stream);
  void clip(T_data min, T_data max) {clip(min, max, 0);};

  /**< transpose */
  int transpose(CarmaObj<T_data> *source);
  // CarmaObj<T_data>& operator= (const CarmaObj<T_data>& obj);

  /*
   *  ____  _        _    ____  _
   * | __ )| |      / \  / ___|/ |
   * |  _ \| |     / _ \ \___ \| |
   * | |_) | |___ / ___ \ ___) | |
   * |____/|_____/_/   \_\____/|_|
   *
   */

  int aimax(int incx);
  int aimin(int incx);
  T_data asum(int incx);
  T_data nrm2(int incx);
  T_data dot(CarmaObj<T_data> *source, int incx, int incy);
  void scale(T_data alpha, int incx);
  void swap(CarmaObj<T_data> *source, int incx, int incy);
  void copy(CarmaObj<T_data> *source, int incx, int incy);
  void axpy(T_data alpha, CarmaObj<T_data> *source, int incx, int incy,
            int offset = 0);
  void rot(CarmaObj<T_data> *source, int incx, int incy, T_data sc, T_data ss);

  /*
   *  ____  _        _    ____ ____
   * | __ )| |      / \  / ___|___ \
   * |  _ \| |     / _ \ \___ \ __) |
   * | |_) | |___ / ___ \ ___) / __/
   * |____/|_____/_/   \_\____/_____|
   *
   */

  void gemv(char trans, T_data alpha, CarmaObj<T_data> *matA, int lda,
            CarmaObj<T_data> *vectx, int incx, T_data beta, int incy);
  void ger(T_data alpha, CarmaObj<T_data> *vectx, int incx,
           CarmaObj<T_data> *vecty, int incy, int lda);
  void symv(char uplo, T_data alpha, CarmaObj<T_data> *matA, int lda,
            CarmaObj<T_data> *vectx, int incx, T_data beta, int incy);

  /*
   *  ____  _        _    ____ _____
   * | __ )| |      / \  / ___|___ /
   * |  _ \| |     / _ \ \___ \ |_ \
   * | |_) | |___ / ___ \ ___) |__) |
   * |____/|_____/_/   \_\____/____/
   *
   */

  void gemm(char transa, char transb, T_data alpha, CarmaObj<T_data> *matA,
            int lda, CarmaObj<T_data> *matB, int ldb, T_data beta, int ldc);
  void symm(char side, char uplo, T_data alpha, CarmaObj<T_data> *matA,
            int lda, CarmaObj<T_data> *matB, int ldb, T_data beta, int ldc);
  void syrk(char uplo, char transa, T_data alpha, CarmaObj<T_data> *matA,
            int lda, T_data beta, int ldc);
  void syrkx(char uplo, char transa, T_data alpha, CarmaObj<T_data> *matA,
             int lda, CarmaObj<T_data> *matB, int ldb, T_data beta, int ldc);
  void geam(char transa, char transb, T_data alpha, CarmaObj<T_data> *matA,
            int lda, T_data beta, CarmaObj<T_data> *matB, int ldb, int ldc);
  void dgmm(char side, CarmaObj<T_data> *matA, int lda,
            CarmaObj<T_data> *vectx, int incx, int ldc);

  /**< Curand */
  int init_prng();
  int init_prng(long seed);
  int destroy_prng();
  int prng(T_data *output, char gtype, float alpha, float beta);
  int prng(T_data *output, char gtype, float alpha);
  int prng(char gtype, float alpha, float beta);
  int prng(char gtype, float alpha);
  int prng(char gtype);

  int prng_montagn(float init_montagn);

  int init_prng_host(int seed);
  int prng_host(char gtype);
  int prng_host(char gtype, T_data stddev);
  int prng_host(char gtype, T_data stddev, T_data alpha);
  int destroy_prng_host();
};
typedef CarmaObj<int> CarmaObjI;
typedef CarmaObj<unsigned int> CarmaObjUI;
typedef CarmaObj<uint16_t> CarmaObjUSI;
typedef CarmaObj<float> CarmaObjS;
typedef CarmaObj<double> CarmaObjD;
typedef CarmaObj<float2> CarmaObjS2;
typedef CarmaObj<double2> CarmaObjD2;
typedef CarmaObj<cuFloatComplex> CarmaObjC;
typedef CarmaObj<cuDoubleComplex> CarmaObjZ;
// typedef CarmaObj<tuple_t<float>> CarmaObjTF;

#ifdef CAN_DO_HALF
typedef CarmaObj<half> CarmaObjH;
#endif

template <class T_data>
std::ostream &operator<<(std::ostream &os, CarmaObj<T_data> &obj) {
  os << "-----------------------" << std::endl;
  os << "CarmaObj<" << typeid(T_data).name() << "> object on GPU"
     << obj.get_device() << std::endl;
  long ndims = obj.get_dims(0);
  os << "ndims = " << ndims << std::endl;
  for (long dim = 0; dim < ndims; dim++) {
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
void clip_array(T_data *d_data, T_data min, T_data max, int N,
                CarmaDevice *device, cudaStream_t stream);

// CU functions sum
template <class T_data>
void reduce(int size, int threads, int blocks, T_data *d_idata,
            T_data *d_odata);
template <class T_data>
T_data reduce(T_data *data, int N);

template <class T_data>
void init_reduceCubCU(T_data *&cub_data, size_t &cub_data_size, T_data *data,
                      T_data *&o_data, int N);
template <class T_data>
void reduceCubCU(T_data *cub_data, size_t cub_data_size, T_data *data,
                 T_data *o_data, int N, cudaStream_t stream=0);

// CU functions transpose
template <class T_data>
int transposeCU(T_data *d_idata, T_data *d_odata, long N1, long N2);

// CU functions generic
template <class T_data>
int launch_generic1d(T_data *d_idata, T_data *d_odata, int N,
                     CarmaDevice *device);
template <class T_data>
int launch_generic2d(T_data *d_odata, T_data *d_idata, int N1, int N2);

// CU functions curand
int carma_prng_init(int *seed, const int nb_threads, const int nb_blocks,
                    curandState *state);
template <class T>
int carma_prng_cu(T *results, const int nb_threads, const int nb_blocks,
                  curandState *state, char gtype, int n, float alpha,
                  float beta);
template <class T>
int carma_curand_montagn(curandState *state, T *d_odata, int N,
                         CarmaDevice *device);

// CU functions fft
template <class T_in, class T_out>
cufftType carma_select_plan();
template <class T_in, class T_out>
void carma_initfft(const long *dims_data, cufftHandle *plan, cufftType type_plan);
template <class T_in, class T_out>
int CarmaFFT(T_in *input, T_out *output, int dir, cufftHandle plan);

// CU functions generic
template <class T_data>
int fillindex(T_data *d_odata, T_data *d_idata, int *indx, int N,
              CarmaDevice *device);
template <class T_data>
int fillvalues(T_data *d_odata, T_data *val, int N, CarmaDevice *device);
template <class T>
int getarray2d(T *d_odata, T *d_idata, int x0, int Ncol, int NC, int N,
               CarmaDevice *device);
template <class T>
int fillarray2d(T *d_odata, T *d_idata, int x0, int Ncol, int NC, int N,
                CarmaDevice *device);
template <class T>
int fillarray2d2(T *d_odata, T *d_idata, int x0, int Ncol, int NC, int N,
                 CarmaDevice *device);
template <class T>
int fill_sym_matrix(char src_uplo, T *d_data, int Ncol, int N,
                    CarmaDevice *device);
template <class T>
int carma_plus(T *d_odata, T elpha, int N, CarmaDevice *device);
template <class T>
int carma_plusai(T *d_odata, T *i_data, int i, int sgn, int N,
                 CarmaDevice *device);

// CU functions fftconv
// int fftconv_unpad(float *d_odata, float *d_idata, int fftW, int dataH,
//                   int dataW, int N, int n, int nim);
// int carma_initfftconv(CarmaObjS *data_in, CarmaObjS *kernel_in, CarmaObjS *padded_data,
//                       CarmaObjC *padded_spectrum, int kernelY, int kernelX);
// // CPP functions fftconv
// int carma_fftconv(CarmaObjS *data_out, CarmaObjS *padded_data,
//                   CarmaObjC *padded_spectrum, int kernelY, int kernelX);

#ifdef CAN_DO_HALF
int custom_half_axpy(half alpha, half *source, int incx, int incy, int N,
                     half *dest, CarmaDevice *device);
#endif

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
int extract(T *d_smallimg, const T *d_fullimg, int fullimg_size, int center_pos,
            int extract_size, bool roll);

#endif  // _CARMA_OBJ_H_
