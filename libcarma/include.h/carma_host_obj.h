/**
 * \class carma_host_obj
 *
 * \ingroup libcarma
 *
 * \brief this class provides wrappers to the generic carma host object
 *
 * \author $Author: dg, as $
 *
 * \version $Revision: 1.0 $
 *
 * \date $Date: 2011/01/28$
 *
 */
#ifndef _CARMA_HOST_OBJ_H_
#define _CARMA_HOST_OBJ_H_

#include <iostream>
#include <carma.h>
#include <carma_utils.h>
#include <carma_streams.h>
#include <carma_context.h>

enum MemAlloc {
  MA_MALLOC, MA_PAGELOCK, MA_ZEROCPY, MA_PORTABLE, MA_WRICOMB, MA_GENEPIN
};

#define MEMORY_ALIGNMENT  4096
#define ALIGN_UP(x,size) ( ((size_t)x+(size-1))&(~(size-1)) )

using namespace std;

template<class T_data>
class carma_obj;

template<class T_data>
class carma_host_obj {

protected:
  T_data *h_data; ///< Input data
  T_data *data_UA; ///< unpadded input dara for generic pinned mem
  long *dims_data; ///< dimensions of the array
  int nb_elem; ///< number of elments in the array
  MemAlloc mallocType; ///< type of host alloc
  carma_streams *streams;

  void init(const long *dims_data,const T_data *data, MemAlloc mallocType,
      int nb_streams);

public:
  carma_host_obj(const long *dims_data);
  carma_host_obj(const long *dims_data, MemAlloc mallocType);
  carma_host_obj(const carma_host_obj<T_data> *obj);
  carma_host_obj(const carma_host_obj<T_data> *obj, MemAlloc mallocType);
  carma_host_obj(const long *dims_data, T_data *data);
  carma_host_obj(const long *dims_data, T_data *data, MemAlloc mallocType);
  carma_host_obj(const long *dims_data, int nb_streams);
  carma_host_obj(const long *dims_data, MemAlloc mallocType, int nb_streams);
  carma_host_obj(const carma_host_obj<T_data> *obj, int nb_streams);
  carma_host_obj(const carma_host_obj<T_data> *obj, MemAlloc mallocType,
      int nb_streams);
  carma_host_obj(const long *dims_data, T_data *data, int nb_streams);
  carma_host_obj(const long *dims_data, T_data *data, MemAlloc mallocType,
      int nb_streams);
  ~carma_host_obj();

  void get_devpntr(void **pntr_dev);

  int get_nbStreams();
  int add_stream();
  int add_stream(int nb);
  int del_stream();
  int del_stream(int nb);
  cudaStream_t get_cudaStream_t(int stream);
  int wait_stream(int stream);
  int wait_all_streams();

  int cpy_obj(carma_obj<T_data>* caObj, cudaMemcpyKind flag);
  int cpy_obj(carma_obj<T_data>* caObj, cudaMemcpyKind flag,
      unsigned int stream);

  /**< General Utilities */
  operator T_data*() {
    return h_data;
  }
  T_data* operator[](int index) {
    return &h_data[index];
  }
  T_data* getData() {
    return h_data;
  }
  T_data* getData(int index) {
    return &h_data[index];
  }
  const long *getDims() {
    return dims_data;
  }
  long getDims(int i) {
    return dims_data[i];
  }
  int getNbElem() {
    return nb_elem;
  }

  /**< Memory transfer */
  int fill_from(const T_data *data);
  int fill_into(T_data *data);

  string getMetAlloc() {
    switch (mallocType) {
    case MA_MALLOC:
      return "MA_MALLOC";
    case MA_PAGELOCK:
      return "MA_PAGELOCK";
    case MA_ZEROCPY:
      return "MA_ZEROCPY";
    case MA_PORTABLE:
      return "MA_PORTABLE";
    case MA_WRICOMB:
      return "MA_WRICOMB";
    case MA_GENEPIN:
      return "MA_GENEPIN";
    default:
      return "MA_UNKNOWN";
    }
  }
};

// MAGMA functions
template<class T_data>
int carma_svd_cpu(carma_host_obj<T_data> *imat,
    carma_host_obj<T_data> *eigenvals, carma_host_obj<T_data> *mod2act,
    carma_host_obj<T_data> *mes2mod);
template<class T>
int carma_getri_cpu(carma_host_obj<T> *h_A);
template<class T>
int carma_potri_cpu(carma_host_obj<T> *h_A);
template<class T>
int carma_syevd_cpu(char jobz, carma_host_obj<T> *h_A,
    carma_host_obj<T> *eigenvals);

// MAGMA functions (direct access)
template<class T>
int carma_svd_cpu(long N, long M, T *imat, T *eigenvals, T *mod2act,
    T *mes2mod);
template<class T>
int carma_getri_cpu(long N, T *h_A);
template<class T>
int carma_potri_cpu(long N, T *h_A);
template<class T>
int carma_syevd_cpu(char jobz, int N, T *h_A, T *eigenvals);
template<class T>
int carma_axpy_cpu(long N, T alpha, T *h_X, long incX, T *h_Y, long incY);
template<class T>
int carma_gemm_cpu(char transa, char transb, long m, long n, long k, T alpha,
    T *A, long lda, T *B, long ldb, T beta, T *C, long ldc);

// CULA functions
template<class T_data>
int carma_cula_svd(carma_host_obj<T_data> *imat,
    carma_host_obj<T_data> *eigenvals, carma_host_obj<T_data> *mod2act,
    carma_host_obj<T_data> *mes2mod);
/* 
 extern "C" {

 }
 */

#endif // _CARMA_HOST_OBJ_H_
