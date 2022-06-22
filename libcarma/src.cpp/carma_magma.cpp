// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      carma_magma.cpp
//! \ingroup   libcarma
//! \brief     this file provides wrappers to the magma functions
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.3.0
//! \date      2022/01/24

#include <carma_magma.h>

#include <type_list.hpp>

#ifdef CAN_DO_HALF
using TypeListObj = GenericTypeList<uint16_t, int, unsigned int, float, double,
                                    half, cuFloatComplex,
                                    cuDoubleComplex>;  // , tuple_t<float>>;
#else
using TypeListObj =
    GenericTypeList<uint16_t, int, unsigned int, float, double, cuFloatComplex,
                    cuDoubleComplex>;  // , tuple_t<float>>;
#endif

#define MAGMA_TRACE(fmt, args...) \
  fprintf(stderr, "%s:%d Warning: " fmt, __FILE__, __LINE__, ##args)

#ifdef USE_MAGMA
#include "magma.h"
#include "magma_lapack.h"
#ifdef USE_MAGMA_SPARSE
#include "magmasparse.h"
#define TEST_USE_SMAGMA(...) __VA_ARGS__
#else  // ifdef USE_MAGMA_SPARSE
// #warning "SPARSE MAGMA will not be used"
#define TEST_USE_SMAGMA(...)                                \
  MAGMA_TRACE("!!!!!! SPARSE MAGMA not compiled !!!!!!\n"); \
  return EXIT_FAILURE;

#endif  // ifdef USE_MAGMA_SPARSE

#if (MAGMA_VERSION_MAJOR == 1) && (MAGMA_VERSION_MINOR == 4)
#define lapacke_vec_const(var) var
#define magma_vec_const(var) var
#define magma_uplo_const(var) var
#endif  // if (MAGMA_VERSION_MAJOR == 1) && (MAGMA_VERSION_MINOR == 4)

#define CHECK_MAGMA(fct, info, dealloc)                                \
  fct;                                                                 \
  if (info) {                                                          \
    printf("%s@%d magma returned error %d: %s.\n", __FILE__, __LINE__, \
           (int)info, magma_strerror(info));                           \
    dealloc;                                                           \
    return EXIT_FAILURE;                                               \
  }

#ifndef max
#define max(a, b) (((a) < (b)) ? (b) : (a))
#endif  // ifndef max
#ifndef min
#define min(a, b) (((a) > (b)) ? (b) : (a))
#endif  // ifndef min

/*
 * Generics template functions
 */

template <class T>
T *create_padded_data(magma_int_t N, const magma_int_t size, magma_int_t &ldda,
                      T *d_src_data) {
  ldda = ((N + (size - 1)) / size) * size;

  //  MAGMA_TRACE(
  //      "the dimension of the matrix (N=%ld) is not a multiple of %ld ->
  //      padding leading dimension
  // to %ld\n",
  //      (long)N, (long)size, (long)ldda);
  T *d_padded_data;
  magma_malloc((void **)&d_padded_data, N * ldda * sizeof(T));
  cudaMemcpy2D(d_padded_data, ldda * sizeof(T), d_src_data, N * sizeof(T),
               N * sizeof(T), N, cudaMemcpyDeviceToDevice);
  return d_padded_data;
}

template <class T>
void destroy_padded_data(magma_int_t N, magma_int_t ldda, T *d_padded_data,
                         T *d_src_data) {
  cudaMemcpy2D(d_src_data, N * sizeof(T), d_padded_data, ldda * sizeof(T),
               N * sizeof(T), N, cudaMemcpyDeviceToDevice);
  magma_free(d_padded_data);
}

template <class T, typename Fn>
int carma_magma_syevd_gen(Fn const &ptr_syevd_gpu, magma_vec_t jobz,
                          magma_int_t N, magma_int_t ldda, T *d_mat,
                          T *h_eigenvals) {
  magma_int_t *iwork;

  //  if (ldda & 31) {
  //    cerr << "Leading dimension of the matrix must be a multiple of 32" <<
  //    std::endl; return EXIT_FAILURE;
  //  }

  magma_int_t info = 0, lwork, liwork, lda = N, aux_iwork[1];
  T aux_work[1];
  T *h_R, *h_work;

  /* Query for workspace sizes */
  CHECK_MAGMA(ptr_syevd_gpu(jobz, magma_uplo_const('L'), N, NULL, ldda, NULL,
                            NULL, lda, aux_work, -1, aux_iwork, -1, &info),
              info, info = 0);

  lwork = (magma_int_t)aux_work[0];
  liwork = aux_iwork[0];
  magma_malloc_cpu((void **)&iwork, (liwork) * sizeof(magma_int_t));

  magma_malloc_pinned((void **)&h_R, (N * lda) * sizeof(T));
  magma_malloc_pinned((void **)&h_work, (lwork) * sizeof(T));

  CHECK_MAGMA(
      ptr_syevd_gpu(jobz, magma_uplo_const('L'), N, d_mat, ldda, h_eigenvals,
                    h_R, lda, h_work, lwork, iwork, liwork, &info),
      info, info = 0);

  magma_free_pinned(h_R);
  magma_free_pinned(h_work);
  magma_free_cpu(iwork);

  return EXIT_SUCCESS;
}

template <class T, typename Fn>
int carma_magma_syevd_m_gen(Fn const &ptr_syevd_m, magma_int_t ngpu,
                            magma_vec_t jobz, magma_int_t N, T *mat,
                            T *eigenvals) {
  magma_int_t *iwork;
  magma_int_t info = 0, lwork, liwork, lda = N, aux_iwork[1];
  T aux_work[1];
  T *h_work;

  /* Query for workspace sizes */
  CHECK_MAGMA(ptr_syevd_m(ngpu, jobz, magma_uplo_const('L'), N, NULL, lda, NULL,
                          aux_work, -1, aux_iwork, -1, &info),
              info, info = 0);

  lwork = (magma_int_t)aux_work[0];
  liwork = aux_iwork[0];

  magma_malloc_cpu((void **)&iwork, (liwork) * sizeof(magma_int_t));

  // magma_malloc_pinned((void**) &h_R, (N * lda) * sizeof(T));
  magma_malloc_pinned((void **)&h_work, (lwork) * sizeof(T));

  CHECK_MAGMA(ptr_syevd_m(ngpu, jobz, magma_uplo_const('L'), N, mat, lda,
                          eigenvals, h_work, lwork, iwork, liwork, &info),
              info, info = 0);

  // magma_free_pinned(h_R);
  magma_free_pinned(h_work);
  magma_free_cpu(iwork);

  return EXIT_SUCCESS;
}

template <class T, typename Fn>
int carma_magma_svd_cpu_gen(Fn const &ptr_gesvd, CarmaHostObj<T> *mat,
                            CarmaHostObj<T> *eigenvals,
                            CarmaHostObj<T> *mod2act,
                            CarmaHostObj<T> *mes2mod) {
  const long *dims = mat->get_dims();
  long m = dims[1];
  long n = dims[2];
  long min_mn = m < n ? m : n;
  long max_mn = m > n ? m : n;

  long *dims_data = new long[2];

  dims_data[0] = 1;
  long nb = 128;  // magma_get_sgesvd_nb(n);
  magma_int_t lwork = max(5 * min_mn, (3 * min_mn + max_mn)) * nb;
  dims_data[1] = lwork;
  CarmaHostObj<T> *h_work =
      new CarmaHostObj<T>(dims_data, MA_PAGELOCK);  // PAGELOCK

  CarmaHostObj<T> *tmp = new CarmaHostObj<T>(mat, MA_PAGELOCK);

  magma_int_t info = 0;
  CHECK_MAGMA(
      ptr_gesvd(magma_vec_const('A'), magma_vec_const('A'), m, n, *tmp, m,
                *eigenvals, *mes2mod, m, *mod2act, n, *h_work, lwork, &info),
      info, delete h_work;
      delete tmp; delete[] dims_data);

  delete h_work;
  delete tmp;
  delete[] dims_data;

  return EXIT_SUCCESS;
}

template <class T, typename Fn>
int carma_magma_potr_inv_gen(Fn const &ptr_potrf, Fn const &ptr_potri,
                          magma_int_t N, magma_int_t ldda, T *d_iA,
                          CarmaDevice *device) {
  magma_int_t info = 0;

  CHECK_MAGMA(ptr_potrf(magma_uplo_const('L'), N, d_iA, ldda, &info), info,
              info = 0);
  CHECK_MAGMA(ptr_potri(magma_uplo_const('L'), N, d_iA, ldda, &info), info,
              info = 0);

  fill_sym_matrix<T>('L', d_iA, N, N * ldda, device);
  return EXIT_SUCCESS;
}

template <class T, typename Fn>
int carma_magma_syevd_cpu_gen(Fn const &ptr_syevd_cpu, magma_vec_t jobz,
                              magma_int_t N, magma_int_t lda, T *h_mat,
                              T *h_eigenvals) {
  magma_int_t *iwork;

  magma_int_t info = 0, liwork = -1, aux_iwork[1];
  T *h_R;

  /* Query for workspace sizes */
  char uplo = 'L';
  char job2 = lapacke_vec_const(jobz);

  CHECK_MAGMA(ptr_syevd_cpu(&job2, &uplo, &N, NULL, &lda, NULL, NULL, &lda,
                            aux_iwork, &liwork, &info),
              info, info = 0);

  liwork = aux_iwork[0];
  magma_malloc_cpu((void **)&iwork, (liwork) * sizeof(magma_int_t));

  cudaMallocHost((void **)&h_R, (N * lda) * sizeof(T));

  CHECK_MAGMA(ptr_syevd_cpu(&job2, &uplo, &N, h_mat, &lda, h_eigenvals, h_R,
                            &lda, iwork, &liwork, &info),
              info, cudaFreeHost(h_R);
              magma_free_cpu(iwork););

  cudaFreeHost(h_R);
  magma_free_cpu(iwork);

  return EXIT_SUCCESS;
}

template <class T, typename Fnf, typename Fni, typename Fni_nb>
int carma_magma_getri_gen(Fnf const &ptr_getrf, Fni const &ptr_getri,
                          Fni_nb const &ptr_get_getri_nb, long N, T *d_iA) {
  magma_int_t *ipiv;

  magma_malloc_cpu((void **)&ipiv, N * sizeof(magma_int_t));

  T *dwork;
  magma_int_t lda = N;
  magma_int_t ldda = ((N + 31) / 32) * 32;
  magma_int_t ldwork = N * ptr_get_getri_nb(N);

  if (lda * ldda != N * N)
    std::cerr << "Please provide arrays of proper size" << std::endl;
  ;

  magma_int_t info = 0;
  cudaMalloc((void **)&dwork, ldwork * sizeof(T));
  CHECK_MAGMA(ptr_getrf(N, N, d_iA, ldda, ipiv, &info), info, cudaFree(dwork);
              magma_free_cpu(ipiv));
  CHECK_MAGMA(ptr_getri(N, d_iA, ldda, ipiv, dwork, ldwork, &info), info,
              cudaFree(dwork);
              magma_free_cpu(ipiv));

  cudaFree(dwork);
  magma_free_cpu(ipiv);
  return EXIT_SUCCESS;
}

template <class T, typename Fnf, typename Fni, typename Fni_nb>
int carma_magma_getri_cpu_gen(Fnf const &ptr_getrf, Fni const &ptr_getri,
                              Fni_nb const &ptr_get_getri_nb, magma_int_t N,
                              T *h_A) {
  magma_int_t *ipiv;

  magma_malloc_cpu((void **)&ipiv, N * sizeof(magma_int_t));

  T *work;
  magma_int_t lda = N;
  magma_int_t lwork = N * ptr_get_getri_nb(N);
  magma_malloc_cpu((void **)&work, lwork * sizeof(T));

  magma_int_t info = 0;
  CHECK_MAGMA(ptr_getrf(&N, &N, h_A, &lda, ipiv, &info), info,
              magma_free_cpu(work);
              magma_free_cpu(ipiv));
  CHECK_MAGMA(ptr_getri(&N, h_A, &lda, ipiv, work, &lwork, &info), info,
              magma_free_cpu(work);
              magma_free_cpu(ipiv));

  magma_free_cpu(work);
  magma_free_cpu(ipiv);
  return EXIT_SUCCESS;
}

template <class T, typename Fnf, typename Fni, typename Fnf_nb,
          typename FnSetMat, typename FnGetMat>
int carma_magma_potr_inv_m_gen(Fnf const &ptr_potrf, Fni const &ptr_potri,
                            Fnf_nb const &ptr_get_potrf_nb,
                            FnSetMat const &ptr_setmatrix_1D_col_bcyclic,
                            FnGetMat const &ptr_getmatrix_1D_col_bcyclic,
                            long num_gpus, T *h_A, T *d_iA, long N,
                            CarmaDevice *device) {
  magma_int_t nb = ptr_get_potrf_nb(N);
  magma_int_t ldda = ((N + nb - 1) / nb) * nb;
  magma_int_t lda = N;

  T *d_lA[num_gpus];
  magma_int_t max_size =
      nb * (1 + N / (nb * num_gpus)) * nb * ((N + nb - 1) / nb);

  for (long dev = 0; dev < num_gpus; dev++) {
    magma_setdevice(dev);
    cudaMalloc(&d_lA[dev], max_size * sizeof(T));
  }
  ptr_setmatrix_1D_col_bcyclic(N, N, h_A, lda, d_lA, ldda, num_gpus, nb);

  magma_int_t info = 0;
  CHECK_MAGMA(
      ptr_potrf(num_gpus, magma_uplo_const('L'), N, d_lA, N, &info), info,

      for (long dev = 0; dev < num_gpus; dev++) {
        magma_setdevice(dev);
        cudaFree(d_lA[dev]);
      });
  CHECK_MAGMA(
      ptr_getmatrix_1D_col_bcyclic(N, N, d_lA, ldda, h_A, lda, num_gpus, nb),
      info,

      for (long dev = 0; dev < num_gpus; dev++) {
        magma_setdevice(dev);
        cudaFree(d_lA[dev]);
      });

  cudaMemcpy(h_A, d_iA, N * N * sizeof(T), cudaMemcpyHostToDevice);

  // d_iA->host2device(h_A);
  CHECK_MAGMA(
      ptr_potri(magma_uplo_const('L'), N, d_iA, N, &info), info,

      for (long dev = 0; dev < num_gpus; dev++) {
        magma_setdevice(dev);
        cudaFree(d_lA[dev]);
      });

  for (long dev = 0; dev < num_gpus; dev++) {
    magma_setdevice(dev);
    cudaFree(d_lA[dev]);
  }

  fill_sym_matrix<T>('L', d_iA, N, N * N, device);
  return EXIT_SUCCESS;
}

template <class T, typename Fn, typename Fn_copy>
int carma_magma_potr_inv_cpu_gen(Fn const &ptr_potrf, Fn const &ptr_potri,
                              Fn_copy const &ptr_copy, magma_int_t N, T *h_A) {
  magma_int_t info = 0;
  char uplo = 'L';
  magma_int_t one = 1;

  CHECK_MAGMA(ptr_potrf(&uplo, &N, h_A, &N, &info), info, info = 0)
  CHECK_MAGMA(ptr_potri(&uplo, &N, h_A, &N, &info), info, info = 0)

  // Copy lower part to the upper
  for (magma_int_t row = 1; row < N; row++) {
    ptr_copy(&row, h_A + row, &N, h_A + row * N, &one);
  }
  return EXIT_SUCCESS;
}

template <class T, typename Fn>
int carma_magma_gemv_gen(Fn mcompute, char trans, int m, int n, T alpha,
                         T *matA, int lda, T *vectx, int incx, T beta, T *vecty,
                         int incy) {
  magma_trans_t trans2 =
      (trans == SOLVER_EIG_MODE_NOVECTOR || trans == 'n') ? MagmaNoTrans : MagmaTrans;
  DEBUG_TRACE("%c %d", trans, trans2);
  mcompute(trans2, m, n, alpha, matA, lda, vectx, incx, beta, vecty, incy);

  return EXIT_SUCCESS;
}

#ifdef USE_MAGMA_SPARSE
template <class T, class var, typename Fn_mtransfert, typename Fn_mconver,
          typename Fn_mfree>
var carma_magma_csr2ell_gen(Fn_mtransfert const &mtransfer,
                            Fn_mconver const &mconvert, Fn_mfree const &mfree,
                            CarmaSparseObj<T> *dA) {
  var A;

  if (dA->format != "CSR") {
    DEBUG_TRACE("carma_fill_magma_sparse_matrix needs a CSR matrix as input");
    return A;
  }

  A.storage_type = Magma_CSR;
  A.memory_location = Magma_DEV;
  A.fill_mode = MagmaFull;
  A.num_rows = dA->dims_data[1];
  A.num_cols = dA->dims_data[2];
  A.nnz = dA->nz_elem;
  A.dval = dA->d_data;
  A.drow = dA->d_rowind;
  A.dcol = dA->d_colind;
  magma_queue_t queue;
  magma_queue_create(/*devices[ opts->device ],*/ &queue);

  var hA, hB;
  mtransfer(A, &hA, Magma_DEV, Magma_CPU, queue);
  mfree(&A, queue);
  mconvert(hA, &hB, Magma_CSR, Magma_ELL, queue);
  mfree(&hA, queue);
  mtransfer(hB, &A, Magma_CPU, Magma_DEV, queue);
  mfree(&hB, queue);
  magma_queue_destroy(queue);

  // DEBUG_TRACE("nnz %d d_data %p d_rowind %p d_colind %p", A.nnz, A.dval,
  // A.drow, A.dcol);
  dA->format = "ELL";
  dA->nz_elem = A.nnz;
  dA->d_data = A.dval;
  dA->d_rowind = A.drow;
  dA->d_colind = A.dcol;

  return A;
}

template <class T, class varM, class varV, typename Fn>
int carma_magma_spmv_gen(Fn spmv, T alpha, varM dA, CarmaObj<T> *dx, T beta,
                         CarmaObj<T> *dy) {
  varV x, y;

  x.memory_location = Magma_DEV;
  x.num_rows = dx->get_nb_elements();
  x.num_cols = 1;
  x.nnz = x.num_rows;
  x.dval = dx->get_data();
  x.major = MagmaColMajor;

  y.memory_location = Magma_DEV;
  y.num_rows = dy->get_nb_elements();
  y.num_cols = 1;
  y.nnz = x.num_rows;
  y.dval = dy->get_data();
  y.major = MagmaColMajor;

  magma_queue_t queue;
  magma_queue_create(/*devices[ opts->device ],*/ &queue);
  spmv(alpha, dA, x, beta, y, queue);
  magma_queue_destroy(queue);

  return EXIT_SUCCESS;
}

#endif  // ifdef USE_MAGMA_SPARSE

#define TEST_USE_MAGMA(...) __VA_ARGS__

#else  // ifdef USE_MAGMA
#warning "MAGMA will not be used"

#define TEST_USE_MAGMA(...)                          \
  return EXIT_FAILURE;
  // MAGMA_TRACE("!!!!!! MAGMA not compiled !!!!!!\n"); \

#warning "SPARSE MAGMA will not be used"
#define TEST_USE_SMAGMA(...)                                \
  return EXIT_FAILURE;
  // MAGMA_TRACE("!!!!!! SPARSE MAGMA not compiled !!!!!!\n"); \

#endif  // ifdef USE_MAGMA

/*
 * API functions
 */
int carma_magma_disabled() { TEST_USE_MAGMA(return EXIT_SUCCESS); }

template <class T>
int carma_magma_syevd(char jobz, long N, T *mat, T *eigenvals) {
  DEBUG_TRACE("Not implemented for this data type");
  return EXIT_FAILURE;
}

template <>
int carma_magma_syevd<float>(char jobz, long N, float *mat, float *eigenvals) {
  TEST_USE_MAGMA(return carma_magma_syevd_gen(
      magma_ssyevd_gpu, magma_vec_const(jobz), N, N, mat, eigenvals));
}

template <>
int carma_magma_syevd<double>(char jobz, long N, double *mat,
                              double *eigenvals) {
  TEST_USE_MAGMA(return carma_magma_syevd_gen(
      magma_dsyevd_gpu, magma_vec_const(jobz), N, N, mat, eigenvals));
}

template <class T>
int carma_magma_syevd(char jobz, CarmaObj<T> *mat,
                      CarmaHostObj<T> *eigenvals) {
  DEBUG_TRACE("Not implemented for this data type");
  return EXIT_FAILURE;
}

template <>
int carma_magma_syevd<float>(char jobz, CarmaObjS *mat,
                             CarmaHostObj<float> *eigenvals) {
  long N = mat->get_dims(1);

  if (N != mat->get_dims(2)) {
    std::cerr << "Matrix must be symmetric" << std::endl;

    return EXIT_FAILURE;
  }
  TEST_USE_MAGMA(
      return carma_magma_syevd_gen(magma_ssyevd_gpu, magma_vec_const(jobz), N,
                                   N, mat->get_data(), eigenvals->get_data()));
}

template <>
int carma_magma_syevd<double>(char jobz, CarmaObjD *mat,
                              CarmaHostObj<double> *eigenvals) {
  long N = mat->get_dims(1);

  if (N != mat->get_dims(2)) {
    std::cerr << "Matrix must be symmetric" << std::endl;

    return EXIT_FAILURE;
  }
  TEST_USE_MAGMA(
      return carma_magma_syevd_gen(magma_dsyevd_gpu, magma_vec_const(jobz), N,
                                   N, mat->get_data(), eigenvals->get_data()));
}

/*
 template<> int carma_magma_syevd<float, 2>(char jobz, CarmaObjS *mat,
 CarmaHostObj<float> *eigenvals) {
 MAGMA_TRACE("carma_magma_syevd: this method not implemented !\n");
 return EXIT_FAILURE;
 }
 template<> int carma_magma_syevd<double, 2>(char jobz, CarmaObjD *mat,
 CarmaHostObj<double> *eigenvals) {
 return carma_magma_syevd<double,
 magma_dsyevd_gpu_magmablas>(magma_vec_const(jobz), mat, eigenvals);
 }
 template<> int carma_magma_syevd<float, 3>(char jobz, CarmaObjS *mat,
 CarmaHostObj<float> *eigenvals) {
 MAGMA_TRACE("carma_magma_syevd: this method not implemented !\n");
 return EXIT_FAILURE;
 }
 template<> int carma_magma_syevd<double, 3>(char jobz, CarmaObjD *mat,
 CarmaHostObj<double> *eigenvals) {
 return carma_magma_syevd<double, magma_dsyevd_gpu_kblas>(magma_vec_const(jobz),
 mat, eigenvals);
 }
 */

template <class T>
int carma_magma_syevd_m(long ngpu, char jobz, long N, T *mat, T *eigenvals) {
  DEBUG_TRACE("Not implemented for this data type");
  return EXIT_FAILURE;
}

template <>
int carma_magma_syevd_m<float>(long ngpu, char jobz, long N, float *mat,
                               float *eigenvals) {
  TEST_USE_MAGMA(
      magma_int_t N_ = N; return carma_magma_syevd_m_gen(
          magma_ssyevd_m, ngpu, magma_vec_const(jobz), N_, mat, eigenvals));
}

template <>
int carma_magma_syevd_m<double>(long ngpu, char jobz, long N, double *mat,
                                double *eigenvals) {
  TEST_USE_MAGMA(
      magma_int_t N_ = N; return carma_magma_syevd_m_gen(
          magma_dsyevd_m, ngpu, magma_vec_const(jobz), N_, mat, eigenvals));
}

template <class T>
int carma_magma_syevd_m(long ngpu, char jobz, long N, CarmaHostObj<T> *mat,
                        CarmaHostObj<T> *eigenvals) {
  DEBUG_TRACE("Not implemented for this data type");
  return EXIT_FAILURE;
}

template <>
int carma_magma_syevd_m<float>(long ngpu, char jobz, CarmaHostObj<float> *mat,
                               CarmaHostObj<float> *eigenvals) {
  long N = mat->get_dims(1);

  if (N != mat->get_dims(2)) {
    std::cerr << "Matrix must be symmetric" << std::endl;

    return EXIT_FAILURE;
  }
  TEST_USE_MAGMA(return carma_magma_syevd_m_gen(
      magma_ssyevd_m, ngpu, magma_vec_const(jobz), N, mat->get_data(),
      eigenvals->get_data()));
}

template <>
int carma_magma_syevd_m<double>(long ngpu, char jobz,
                                CarmaHostObj<double> *mat,
                                CarmaHostObj<double> *eigenvals) {
  long N = mat->get_dims(1);

  if (N != mat->get_dims(2)) {
    std::cerr << "Matrix must be symmetric" << std::endl;

    return EXIT_FAILURE;
  }
  TEST_USE_MAGMA(return carma_magma_syevd_m_gen(
      magma_dsyevd_m, magma_int_t(ngpu), magma_vec_const(jobz), N,
      mat->get_data(), eigenvals->get_data()));
}

template <class T>
int carma_magma_syevd_m(long ngpu, char jobz, CarmaHostObj<T> *mat,
                        CarmaHostObj<T> *eigenvals, CarmaHostObj<T> *U) {
  DEBUG_TRACE("Not implemented for this data type");
  return EXIT_FAILURE;
}

template <>
int carma_magma_syevd_m<float>(long ngpu, char jobz, CarmaHostObj<float> *mat,
                               CarmaHostObj<float> *eigenvals,
                               CarmaHostObj<float> *U) {
  TEST_USE_MAGMA(

      /* Initialize the matrix */
      magma_int_t N = mat->get_dims(1);

      if (N != mat->get_dims(2)) {
        std::cerr << "Matrix must be symmetric" << std::endl;
        return EXIT_FAILURE;
      }

      magma_int_t lda = N;
      lapackf77_slacpy(MagmaUpperLowerStr, &N, &N, *mat, &lda, *U, &lda);
      return carma_magma_syevd_m_gen<float>(
          magma_ssyevd_m, ngpu, magma_vec_const(jobz), N, *U, *eigenvals));
}

template <>
int carma_magma_syevd_m<double>(long ngpu, char jobz,
                                CarmaHostObj<double> *mat,
                                CarmaHostObj<double> *eigenvals,
                                CarmaHostObj<double> *U) {
  TEST_USE_MAGMA(

      /* Initialize the matrix */
      magma_int_t N = mat->get_dims(1);

      if (N != mat->get_dims(2)) {
        std::cerr << "Matrix must be symmetric" << std::endl;
        return EXIT_FAILURE;
      } magma_int_t lda = N;
      lapackf77_dlacpy(MagmaUpperLowerStr, &N, &N, *mat, &lda, *U, &lda);
      return carma_magma_syevd_m_gen(magma_dsyevd_m, ngpu,
                                     magma_vec_const(jobz), N, U->get_data(),
                                     eigenvals->get_data()));
}

template <class T>
int carma_magma_svd_cpu(CarmaHostObj<T> *mat, CarmaHostObj<T> *eigenvals,
                        CarmaHostObj<T> *mod2act,
                        CarmaHostObj<T> *mes2mod) {
  DEBUG_TRACE("Not implemented for this data type");
  return EXIT_FAILURE;
}

template <>
int carma_magma_svd_cpu<float>(CarmaHostObj<float> *mat,
                               CarmaHostObj<float> *eigenvals,
                               CarmaHostObj<float> *mod2act,
                               CarmaHostObj<float> *mes2mod) {
  TEST_USE_MAGMA(return carma_magma_svd_cpu_gen<float>(
      magma_sgesvd, mat, eigenvals, mod2act, mes2mod));
}

template <>
int carma_magma_svd_cpu<double>(CarmaHostObj<double> *mat,
                                CarmaHostObj<double> *eigenvals,
                                CarmaHostObj<double> *mod2act,
                                CarmaHostObj<double> *mes2mod) {
  TEST_USE_MAGMA(return carma_magma_svd_cpu_gen<double>(
      magma_dgesvd, mat, eigenvals, mod2act, mes2mod));
}

template <class T>
int carma_magma_potr_inv(CarmaObj<T> *d_iA) {
  DEBUG_TRACE("Not implemented for this data type");
  return EXIT_FAILURE;
}

template <>
int carma_magma_potr_inv<float>(CarmaObj<float> *d_iA) {
  const long *dims = d_iA->get_dims();
  long N = dims[1];

  if (N != dims[2]) {
    MAGMA_TRACE(
        "carma_magma_potr_inv : non square and positive-definite matrix\n");

    return EXIT_FAILURE;
  }
  TEST_USE_MAGMA(return carma_magma_potr_inv_gen<float>(
      magma_spotrf_gpu, magma_spotri_gpu, N, N, *d_iA,
      d_iA->get_context()->get_device(d_iA->get_device())));
}

template <>
int carma_magma_potr_inv<double>(CarmaObj<double> *d_iA) {
  const long *dims = d_iA->get_dims();
  long N = dims[1];

  if (N != dims[2]) {
    MAGMA_TRACE(
        "carma_magma_potr_inv : non square and positive-definite matrix\n");

    return EXIT_FAILURE;
  }
  TEST_USE_MAGMA(return carma_magma_potr_inv_gen<double>(
      magma_dpotrf_gpu, magma_dpotri_gpu, N, N, *d_iA,
      d_iA->get_context()->get_device(d_iA->get_device())));
}

template <class T>
int carma_magma_potr_inv_m(long num_gpus, CarmaHostObj<T> *h_A,
                        CarmaObj<T> *d_iA) {
  DEBUG_TRACE("Not implemented for this data type");
  return EXIT_FAILURE;
}

template <>
int carma_magma_potr_inv_m<float>(long num_gpus, CarmaHostObj<float> *h_A,
                               CarmaObj<float> *d_iA) {
  long N = h_A->get_dims(1);

  if (N != h_A->get_dims(2)) {
    std::cerr << "Matrix must be square and positive-definite" << std::endl;

    return EXIT_FAILURE;
  }
  TEST_USE_MAGMA(return carma_magma_potr_inv_m_gen<float>(
      magma_spotrf_mgpu, magma_spotri_gpu, magma_get_spotrf_nb,
      magma_ssetmatrix_1D_col_bcyclic, magma_sgetmatrix_1D_col_bcyclic,
      num_gpus, *h_A, *d_iA, N,
      d_iA->get_context()->get_device(d_iA->get_device())));
}

template <>
int carma_magma_potr_inv_m<double>(long num_gpus, CarmaHostObj<double> *h_A,
                                CarmaObj<double> *d_iA) {
  long N = h_A->get_dims(1);

  if (N != h_A->get_dims(2)) {
    std::cerr << "Matrix must be square and positive-definite" << std::endl;

    return EXIT_FAILURE;
  }

  TEST_USE_MAGMA(return carma_magma_potr_inv_m_gen<double>(
      magma_dpotrf_mgpu, magma_dpotri_gpu, magma_get_dpotrf_nb,
      magma_dsetmatrix_1D_col_bcyclic, magma_dgetmatrix_1D_col_bcyclic,
      num_gpus, *h_A, *d_iA, N,
      d_iA->get_context()->get_device(d_iA->get_device())));
}

template <class T>
int carma_magma_getri(CarmaObj<T> *d_iA) {
  DEBUG_TRACE("Not implemented for this data type");
  return EXIT_FAILURE;
}

template <>
int carma_magma_getri<float>(CarmaObj<float> *d_iA) {
  long N = d_iA->get_dims(1);

  if (N != d_iA->get_dims(2)) {
    std::cerr << "Matrix must be square" << std::endl;

    return EXIT_FAILURE;
  }

  TEST_USE_MAGMA(return carma_magma_getri_gen<float>(
      magma_sgetrf_gpu, magma_sgetri_gpu, magma_get_sgetri_nb, N, *d_iA));
}

template <>
int carma_magma_getri<double>(CarmaObj<double> *d_iA) {
  long N = d_iA->get_dims(1);

  if (N != d_iA->get_dims(2)) {
    std::cerr << "Matrix must be square" << std::endl;

    return EXIT_FAILURE;
  }
  TEST_USE_MAGMA(return carma_magma_getri_gen<double>(
      magma_dgetrf_gpu, magma_dgetri_gpu, magma_get_dgetri_nb, N, *d_iA));
}

template <class T>
int carma_magma_potr_inv_cpu(CarmaHostObj<T> *h_A) {
  DEBUG_TRACE("Not implemented for this data type");
  return EXIT_FAILURE;
}

template <>
int carma_magma_potr_inv_cpu<float>(CarmaHostObj<float> *h_A) {
  long N = h_A->get_dims(1);

  if (N != h_A->get_dims(2)) {
    std::cerr << "Matrix must be symmetric" << std::endl;

    return EXIT_FAILURE;
  }
  TEST_USE_MAGMA(return carma_magma_potr_inv_cpu_gen<float>(
      lapackf77_spotrf, lapackf77_spotri, blasf77_scopy, N, *h_A));
}

template <>
int carma_magma_potr_inv_cpu<double>(CarmaHostObj<double> *h_A) {
  long N = h_A->get_dims(1);

  if (N != h_A->get_dims(2)) {
    std::cerr << "Matrix must be symmetric" << std::endl;

    return EXIT_FAILURE;
  }
  TEST_USE_MAGMA(return carma_magma_potr_inv_cpu_gen<double>(
      lapackf77_dpotrf, lapackf77_dpotri, blasf77_dcopy, N, *h_A));
}

template <class T>
int carma_magma_potr_inv_cpu(long N, T *h_A) {
  DEBUG_TRACE("Not implemented for this data type");
  return EXIT_FAILURE;
}

template <>
int carma_magma_potr_inv_cpu<float>(long N, float *h_A) {
  TEST_USE_MAGMA(return carma_magma_potr_inv_cpu_gen<float>(
      lapackf77_spotrf, lapackf77_spotri, blasf77_scopy, N, h_A));
}

template <>
int carma_magma_potr_inv_cpu<double>(long N, double *h_A) {
  TEST_USE_MAGMA(return carma_magma_potr_inv_cpu_gen<double>(
      lapackf77_dpotrf, lapackf77_dpotri, blasf77_dcopy, N, h_A));
}

template <class T>
int carma_magma_getri_cpu(CarmaHostObj<T> *h_A) {
  DEBUG_TRACE("Not implemented for this data type");
  return EXIT_FAILURE;
}

template <>
int carma_magma_getri_cpu<float>(CarmaHostObj<float> *h_A) {
  long N = h_A->get_dims(1);

  if (N != h_A->get_dims(2)) {
    std::cerr << "Matrix must be symmetric" << std::endl;

    return EXIT_FAILURE;
  }

  TEST_USE_MAGMA(return carma_magma_getri_cpu_gen<float>(
      lapackf77_sgetrf, lapackf77_sgetri, magma_get_dgetri_nb, N, *h_A));
}

template <>
int carma_magma_getri_cpu<double>(CarmaHostObj<double> *h_A) {
  long N = h_A->get_dims(1);

  if (N != h_A->get_dims(2)) {
    std::cerr << "Matrix must be symmetric" << std::endl;

    return EXIT_FAILURE;
  }
  TEST_USE_MAGMA(return carma_magma_getri_cpu_gen<double>(
      lapackf77_dgetrf, lapackf77_dgetri, magma_get_sgetri_nb, N, *h_A));
}

template <class T>
int carma_magma_getri_cpu(long N, T *h_A) {
  DEBUG_TRACE("Not implemented for this data type");
  return EXIT_FAILURE;
}

template <>
int carma_magma_getri_cpu<float>(long N, float *h_A) {
  TEST_USE_MAGMA(return carma_magma_getri_cpu_gen<float>(
      lapackf77_sgetrf, lapackf77_sgetri, magma_get_dgetri_nb, N, h_A));
}

template <>
int carma_magma_getri_cpu<double>(long N, double *h_A) {
  TEST_USE_MAGMA(return carma_magma_getri_cpu_gen<double>(
      lapackf77_dgetrf, lapackf77_dgetri, magma_get_sgetri_nb, N, h_A));
}

template <class T>
int carma_magma_syevd_cpu(char jobz, CarmaHostObj<T> *h_A,
                          CarmaHostObj<T> *eigenvals) {
  DEBUG_TRACE("Not implemented for this data type");
  return EXIT_FAILURE;
}

template <>
int carma_magma_syevd_cpu<float>(char jobz, CarmaHostObj<float> *h_A,
                                 CarmaHostObj<float> *eigenvals) {
  long N = h_A->get_dims(1);

  if (N != h_A->get_dims(2)) {
    std::cerr << "Matrix must be symmetric" << std::endl;

    return EXIT_FAILURE;
  }
  TEST_USE_MAGMA(return carma_magma_syevd_cpu_gen<float>(
      lapackf77_ssyevd, magma_vec_const(jobz), N, N, *h_A, *eigenvals));
}

template <>
int carma_magma_syevd_cpu<double>(char jobz, CarmaHostObj<double> *h_A,
                                  CarmaHostObj<double> *eigenvals) {
  long N = h_A->get_dims(1);

  if (N != h_A->get_dims(2)) {
    std::cerr << "Matrix must be symmetric" << std::endl;

    return EXIT_FAILURE;
  }
  TEST_USE_MAGMA(return carma_magma_syevd_cpu_gen<double>(
      lapackf77_dsyevd, magma_vec_const(jobz), N, N, *h_A, *eigenvals));
}

template <class T>
int carma_magma_syevd_cpu(char jobz, long N, T *h_A, T *eigenvals) {
  DEBUG_TRACE("Not implemented for this data type");
  return EXIT_FAILURE;
}

template <>
int carma_magma_syevd_cpu<float>(char jobz, long N, float *h_A,
                                 float *eigenvals) {
  TEST_USE_MAGMA(
      magma_int_t N_ = N; return carma_magma_syevd_cpu_gen<float>(
          lapackf77_ssyevd, magma_vec_const(jobz), N_, N_, h_A, eigenvals));
}

template <>
int carma_magma_syevd_cpu<double>(char jobz, long N, double *h_A,
                                  double *eigenvals) {
  TEST_USE_MAGMA(
      magma_int_t N_ = N; return carma_magma_syevd_cpu_gen<double>(
          lapackf77_dsyevd, magma_vec_const(jobz), N_, N_, h_A, eigenvals));
}

template <class T>
int carma_magma_axpy_cpu(long N, T alpha, T *h_X, long incX, T *h_Y,
                         long incY) {
  DEBUG_TRACE("Not implemented for this data type");
  return EXIT_FAILURE;
}

template <>
int carma_magma_axpy_cpu<float>(long N, float alpha, float *h_X, long incX,
                                float *h_Y, long incY) {
  TEST_USE_MAGMA(magma_int_t tmp_N = N; magma_int_t tmp_incX = incX;
                 magma_int_t tmp_incY = incY;
                 blasf77_saxpy(&tmp_N, &alpha, h_X, &tmp_incX, h_Y, &tmp_incY));
  return EXIT_SUCCESS;
}

template <>
int carma_magma_axpy_cpu<double>(long N, double alpha, double *h_X, long incX,
                                 double *h_Y, long incY) {
  TEST_USE_MAGMA(magma_int_t tmp_N = N; magma_int_t tmp_incX = incX;
                 magma_int_t tmp_incY = incY;
                 blasf77_daxpy(&tmp_N, &alpha, h_X, &tmp_incX, h_Y, &tmp_incY));

  return EXIT_SUCCESS;
}

template <class T>
int carma_gemm_cpu(char transa, char transb, long m, long n, long k, T alpha,
                   T *A, long lda, T *B, long ldb, T beta, T *C, long ldc) {
  DEBUG_TRACE("Not implemented for this data type");
  return EXIT_FAILURE;
}

template <>
int carma_gemm_cpu<float>(char transa, char transb, long m, long n, long k,
                          float alpha, float *A, long lda, float *B, long ldb,
                          float beta, float *C, long ldc) {
  TEST_USE_MAGMA(magma_int_t tmp_m = m; magma_int_t tmp_n = n;
                 magma_int_t tmp_k = k; magma_int_t tmp_lda = lda;
                 magma_int_t tmp_ldb = ldb; magma_int_t tmp_ldc = ldc;
                 blasf77_sgemm(&transa, &transb, &tmp_m, &tmp_n, &tmp_k, &alpha,
                               A, &tmp_lda, B, &tmp_ldb, &beta, C, &tmp_ldc));
  return EXIT_SUCCESS;
}

template <>
int carma_gemm_cpu<double>(char transa, char transb, long m, long n, long k,
                           double alpha, double *A, long lda, double *B,
                           long ldb, double beta, double *C, long ldc) {
  TEST_USE_MAGMA(magma_int_t tmp_m = m; magma_int_t tmp_n = n;
                 magma_int_t tmp_k = k; magma_int_t tmp_lda = lda;
                 magma_int_t tmp_ldb = ldb; magma_int_t tmp_ldc = ldc;
                 blasf77_dgemm(&transa, &transb, &tmp_m, &tmp_n, &tmp_k, &alpha,
                               A, &tmp_lda, B, &tmp_ldb, &beta, C, &tmp_ldc));
  return EXIT_SUCCESS;
}

template <class T>
int carma_magma_gemv(char trans, int m, int n, T alpha, T *matA, int lda,
                     T *vectx, int incx, T beta, T *vecty, int incy) {
  DEBUG_TRACE("Not implemented for this data type");
  return EXIT_FAILURE;
}

template <>
int carma_magma_gemv<float>(char trans, int m, int n, float alpha, float *matA,
                            int lda, float *vectx, int incx, float beta,
                            float *vecty, int incy) {
  TEST_USE_MAGMA(
      return carma_magma_gemv_gen<float>(magma_sgemv, trans, m, n, alpha, matA,
                                         lda, vectx, incx, beta, vecty, incy);

  );
}

template <>
int carma_magma_gemv<double>(char trans, int m, int n, double alpha,
                             double *matA, int lda, double *vectx, int incx,
                             double beta, double *vecty, int incy) {
  TEST_USE_MAGMA(
      return carma_magma_gemv_gen<double>(magma_dgemv, trans, m, n, alpha, matA,
                                          lda, vectx, incx, beta, vecty, incy);

  );
}

template <class T>
int carma_magma_csr2ell(CarmaSparseObj<T> *dA) {
  DEBUG_TRACE("Not implemented for this data type");
  return EXIT_FAILURE;
}

template <>
int carma_magma_csr2ell<float>(CarmaSparseObj<float> *dA) {
  TEST_USE_SMAGMA(
      dA->s_sparse_mat = carma_magma_csr2ell_gen<float, magma_s_sparse_matrix>(
          magma_s_mtransfer, magma_s_mconvert, magma_s_mfree, dA);

      return EXIT_SUCCESS;

  );
}

template <>
int carma_magma_csr2ell<double>(CarmaSparseObj<double> *dA) {
  TEST_USE_SMAGMA(
      dA->d_sparse_mat = carma_magma_csr2ell_gen<double, magma_d_sparse_matrix>(
          magma_d_mtransfer, magma_d_mconvert, magma_d_mfree, dA);
      return EXIT_SUCCESS;

  );
}

template <class T>
int carma_magma_spmv(T alpha, CarmaSparseObj<T> *dA, CarmaObj<T> *dx, T beta,
                     CarmaObj<T> *dy) {
  DEBUG_TRACE("Not implemented for this data type");
  return EXIT_FAILURE;
}

template <>
int carma_magma_spmv<float>(float alpha, CarmaSparseObj<float> *dA,
                            CarmaObj<float> *dx, float beta,
                            CarmaObj<float> *dy) {
  if (dA->format != "ELL") {
    DEBUG_TRACE("carma_fill_magma_sparse_matrix needs a ELL matrix as input");
    return EXIT_FAILURE;
  }
  TEST_USE_SMAGMA(
      return carma_magma_spmv_gen<float, magma_s_sparse_matrix, magma_s_vector>(
                 magma_s_spmv, alpha, dA->s_sparse_mat, dx, beta, dy);

  );
}

template <>
int carma_magma_spmv<double>(double alpha, CarmaSparseObj<double> *dA,
                             CarmaObj<double> *dx, double beta,
                             CarmaObj<double> *dy) {
  if (dA->format != "ELL") {
    DEBUG_TRACE("carma_fill_magma_sparse_matrix needs a ELL matrix as input");
    return EXIT_FAILURE;
  }
  TEST_USE_SMAGMA(return carma_magma_spmv_gen<double, magma_d_sparse_matrix,
                                              magma_d_vector>(
                             magma_d_spmv, alpha, dA->d_sparse_mat, dx, beta, dy);

  );
}

template <class T>
int carma_magma_sparse_free(CarmaSparseObj<T> *dA) {
  DEBUG_TRACE("Not implemented for this data type");
  return EXIT_FAILURE;
}

template <>
int carma_magma_sparse_free<float>(CarmaSparseObj<float> *dA) {
  if (dA->format != "ELL") {
    return EXIT_FAILURE;
  }
  TEST_USE_SMAGMA(magma_queue_t queue;
                  magma_queue_create(/*devices[ opts->device ],*/ &queue);
                  int ret = magma_s_mfree(&(dA->s_sparse_mat), queue);
                  magma_queue_destroy(queue);

                  dA->format = "NONE"; dA->nz_elem = dA->s_sparse_mat.nnz;
                  dA->d_data = dA->s_sparse_mat.dval;
                  dA->d_rowind = dA->s_sparse_mat.drow;
                  dA->d_colind = dA->s_sparse_mat.dcol; return ret;

  );
}

template <>
int carma_magma_sparse_free<double>(CarmaSparseObj<double> *dA) {
  if (dA->format != "ELL") {
    return EXIT_FAILURE;
  }
  TEST_USE_SMAGMA(magma_queue_t queue;
                  magma_queue_create(/*devices[ opts->device ],*/ &queue);
                  int ret = magma_d_mfree(&(dA->d_sparse_mat), queue);
                  magma_queue_destroy(queue);

                  dA->format = "NONE"; dA->nz_elem = dA->d_sparse_mat.nnz;
                  dA->d_data = dA->d_sparse_mat.dval;
                  dA->d_rowind = dA->d_sparse_mat.drow;
                  dA->d_colind = dA->d_sparse_mat.dcol; return ret;

  );
}

struct CarmaMagmaInterfacer {
  template <typename T_data>
  static void call() {
    force_keep((int (*)(char, long, T_data *, T_data *)) &
               carma_magma_syevd<T_data>);
    force_keep((int (*)(char, CarmaObj<T_data> *, CarmaHostObj<T_data> *)) &
               carma_magma_syevd<T_data>);
    force_keep((int (*)(long, char, long, T_data *, T_data *)) &
               carma_magma_syevd_m<T_data>);
    force_keep((int (*)(long, char, long, CarmaHostObj<T_data> *,
                        CarmaHostObj<T_data> *)) &
               carma_magma_syevd_m<T_data>);
    force_keep((int (*)(long, char, CarmaHostObj<T_data> *,
                        CarmaHostObj<T_data> *, CarmaHostObj<T_data> *)) &
               carma_magma_syevd_m<T_data>);
    force_keep(&carma_magma_svd_cpu<T_data>);
    force_keep(&carma_magma_getri<T_data>);
    force_keep(&carma_magma_potr_inv<T_data>);
    force_keep(&carma_magma_potr_inv_m<T_data>);
    force_keep((int (*)(CarmaHostObj<T_data> *)) &
               carma_magma_potr_inv_cpu<T_data>);
    force_keep((int (*)(long, T_data *)) & carma_magma_potr_inv_cpu<T_data>);
    force_keep((int (*)(CarmaHostObj<T_data> *)) &
               carma_magma_getri_cpu<T_data>);
    force_keep((int (*)(long, T_data *)) & carma_magma_getri_cpu<T_data>);
    force_keep(
        (int (*)(char, CarmaHostObj<T_data> *, CarmaHostObj<T_data> *)) &
        carma_magma_syevd_cpu<T_data>);
    force_keep((int (*)(char, long, T_data *, T_data *)) &
               carma_magma_syevd_cpu<T_data>);
    force_keep(&carma_magma_gemv<T_data>);
    force_keep(&carma_magma_csr2ell<T_data>);
    force_keep(&carma_magma_spmv<T_data>);
    force_keep(&carma_magma_sparse_free<T_data>);
  }
};

void declare_carma_magma() { apply<CarmaMagmaInterfacer, TypeListObj>(); }
