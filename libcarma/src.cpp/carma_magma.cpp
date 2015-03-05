#include <carma_obj.h>
#include <carma_host_obj.h>

#define MAGMA_TRACE(fmt, args...) fprintf(stderr, "%s:%d Warning: " fmt, __FILE__, __LINE__, ## args)

#ifdef USE_MAGMA
#include "magma.h"
#include "magma_lapack.h"

#if (MAGMA_VERSION_MAJOR == 1) && (MAGMA_VERSION_MINOR == 4)
#define lapacke_vec_const(var) var
#define magma_vec_const(var) var
#define magma_uplo_const(var) var
#endif

#define CHECK_MAGMA(fct, info) fct; \
    if (info != 0){							\
      printf("%s@%d magma returned error %d: %s.\n", __FILE__, __LINE__,	\
       (int) info, magma_strerror(info));				\
      return EXIT_FAILURE;						\
    }

#ifndef max
#define max(a,b)  (((a)<(b))?(b):(a))
#endif
#ifndef min
#define min(a,b)  (((a)>(b))?(b):(a))
#endif

/*
 * Generics template functions
 */

template<class T>
T* create_padded_data(magma_int_t N, const magma_int_t size, magma_int_t &ldda,
                      T* d_src_data) {
  ldda = ((N + (size - 1)) / size) * size;
  //  MAGMA_TRACE(
  //      "the dimension of the matrix (N=%ld) is not a multiple of %ld -> padding leading dimension to %ld\n",
  //      (long)N, (long)size, (long)ldda);
  T *d_padded_data;
  magma_malloc((void**) &d_padded_data, N * ldda * sizeof(T));
  cudaMemcpy2D(d_padded_data, ldda * sizeof(T), d_src_data, N * sizeof(T),
               N * sizeof(T), N, cudaMemcpyDeviceToDevice);
  return d_padded_data;
}
template<class T>
void destroy_padded_data(magma_int_t N, magma_int_t ldda, T* d_padded_data,
                         T* d_src_data) {
  cudaMemcpy2D(d_src_data, N * sizeof(T), d_padded_data, ldda * sizeof(T),
               N * sizeof(T), N, cudaMemcpyDeviceToDevice);
  magma_free(d_padded_data);
}

template<class T, magma_int_t (*ptr_syevd_gpu)(magma_vec_t jobz,
                                               magma_uplo_t uplo, magma_int_t n,
                                               T *da, magma_int_t ldda, T *w,
                                               T *wa, magma_int_t ldwa, T *work,
                                               magma_int_t lwork,
                                               magma_int_t *iwork,
                                               magma_int_t liwork,
                                               magma_int_t *info)>
int carma_syevd_gen(magma_vec_t jobz, magma_int_t N, magma_int_t ldda, T *d_mat,
                    T *h_eigenvals) {
  magma_int_t *iwork;
  //  if (ldda & 31) {
  //    cerr << "Leading dimension of the matrix must be a multiple of 32\n";
  //    return EXIT_FAILURE;
  //  }

  magma_int_t info = 0, lwork, liwork, lda = N, aux_iwork[1];
  T aux_work[1];
  T *h_R, *h_work;

  /* Query for workspace sizes */
  CHECK_MAGMA(
      ptr_syevd_gpu(jobz, magma_uplo_const('L'), N, NULL, ldda, NULL, NULL, lda, aux_work, -1, aux_iwork, -1, &info),
      info);

  lwork = (magma_int_t) aux_work[0];
  liwork = aux_iwork[0];
  magma_malloc_cpu((void**) &iwork, (liwork) * sizeof(magma_int_t));

  cudaMallocHost((void**) &h_R, (N * lda) * sizeof(T));
  cudaMallocHost((void**) &h_work, (lwork) * sizeof(T));

  CHECK_MAGMA(
      ptr_syevd_gpu(jobz, magma_uplo_const('L'), N, d_mat, ldda, h_eigenvals,
                    h_R, lda, h_work, lwork, iwork, liwork, &info),
      info);

  cudaFreeHost(h_R);
  cudaFreeHost(h_work);
  magma_free_cpu(iwork);

  return EXIT_SUCCESS;
}

template<class T, magma_int_t (*ptr_syevd_m)(magma_int_t nrgpu,
                                             magma_vec_t jobz,
                                             magma_uplo_t uplo, magma_int_t n,
                                             T *a, magma_int_t lda, T *w,
                                             T *work, magma_int_t lwork,
                                             magma_int_t *iwork,
                                             magma_int_t liwork,
                                             magma_int_t *info)>
int carma_syevd_m_gen(magma_int_t ngpu, magma_uplo_t jobz, magma_int_t N,
                      T *mat, T *eigenvals) {
  magma_int_t *iwork;
  magma_int_t info = 0, lwork, liwork, lda = N, aux_iwork[1];
  T aux_work[1];
  T *h_work;

  /* Query for workspace sizes */
  CHECK_MAGMA(
      ptr_syevd_m(ngpu, magma_vec_const(jobz), magma_uplo_const('L'), N, NULL, lda, NULL, aux_work, -1, aux_iwork, -1, &info),
      info);

  //fprintf(stderr, "%s@%d : I'm here\n", __FILE__, __LINE__);
  lwork = (magma_int_t) aux_work[0];
  liwork = aux_iwork[0];

  magma_malloc_cpu((void**) &iwork, (liwork) * sizeof(magma_int_t));
  //cudaMallocHost((void**) &h_R, (N * lda) * sizeof(T));
  cudaMallocHost((void**) &h_work, (lwork) * sizeof(T));

  CHECK_MAGMA(
      ptr_syevd_m(ngpu, magma_vec_const(jobz), magma_uplo_const('L'), N, mat,
                  lda, eigenvals, h_work, lwork, iwork, liwork, &info),
      info);

  //cudaFreeHost(h_R);
  cudaFreeHost(h_work);
  magma_free_cpu(iwork);

  return EXIT_SUCCESS;
}

template<class T, magma_int_t (*ptr_gesvd)(magma_vec_t jobu, magma_vec_t jobvt,
                                           magma_int_t m, magma_int_t n, T *A,
                                           magma_int_t lda, T *s, T *U,
                                           magma_int_t ldu, T *VT,
                                           magma_int_t ldvt, T *work,
                                           magma_int_t lwork, magma_int_t *info)>
int carma_svd_cpu_gen(carma_host_obj<T> *mat, carma_host_obj<T> *eigenvals,
                      carma_host_obj<T> *mod2act, carma_host_obj<T> *mes2mod) {
  const long *dims = mat->getDims();
  long m = dims[1];
  long n = dims[2];
  long min_mn = m < n ? m : n;
  long max_mn = m > n ? m : n;

  long *dims_data = new long[2];
  dims_data[0] = 1;
  long nb = 128; //magma_get_sgesvd_nb(n);
  magma_int_t lwork = max(5*min_mn, (3*min_mn + max_mn)) * nb;
  dims_data[1] = lwork;
  carma_host_obj<T> *h_work = new carma_host_obj<T>(dims_data, MA_PAGELOCK); //PAGELOCK

  carma_host_obj<T> *tmp = new carma_host_obj<T>(mat, MA_PAGELOCK);

  magma_int_t info = 0;
  CHECK_MAGMA(
      ptr_gesvd(magma_vec_const('A'), magma_vec_const('A'), m, n, *tmp, m,
                *eigenvals, *mes2mod, m, *mod2act, n, *h_work, lwork, &info),
      info);

  delete h_work;
  delete tmp;

  return EXIT_SUCCESS;
}

template<class T, magma_int_t (*ptr_potrf)(magma_uplo_t uplo, magma_int_t n,
                                           T *d_A, magma_int_t ldda,
                                           magma_int_t *info),
    magma_int_t (*ptr_potri)(magma_uplo_t uplo, magma_int_t n, T *d_A,
                             magma_int_t ldda, magma_int_t *info)>
int carma_potri_gen(magma_int_t N, magma_int_t ldda, T *d_iA,
                    carma_device *device) {

  magma_int_t info = 0;
  CHECK_MAGMA(ptr_potrf(magma_uplo_const('L'), N, d_iA, ldda, &info), info);
  CHECK_MAGMA(ptr_potri(magma_uplo_const('L'), N, d_iA, ldda, &info), info);

  fill_sym_matrix<T>('L', d_iA, N, N * ldda, device);
  return EXIT_SUCCESS;
}

template<class T, void (*ptr_syevd_cpu)(const char *jobz, const char *uplo,
                                        const magma_int_t *n, T *da,
                                        const magma_int_t *lda, T *w, T *wa,
                                        const magma_int_t *ldwa,
                                        magma_int_t *iwork,
                                        const magma_int_t *liwork,
                                        magma_int_t *info)>
int carma_syevd_cpu_gen(magma_vec_t jobz, magma_int_t N, magma_int_t lda,
                        T *d_mat, T *h_eigenvals) {
  magma_int_t *iwork;

  magma_int_t info = 0, liwork = -1, aux_iwork[1];
  T *h_R;

  /* Query for workspace sizes */
  char uplo = 'L';
  char job2 = lapacke_vec_const(jobz);
  CHECK_MAGMA(
      ptr_syevd_cpu(&job2, &uplo, &N, NULL, &lda, NULL, NULL, &lda, aux_iwork, &liwork, &info),
      info);

  liwork = aux_iwork[0];
  magma_malloc_cpu((void**) &iwork, (liwork) * sizeof(magma_int_t));

  cudaMallocHost((void**) &h_R, (N * lda) * sizeof(T));

  CHECK_MAGMA(
      ptr_syevd_cpu(&job2, &uplo, &N, d_mat, &lda, h_eigenvals, h_R, &lda,
                    iwork, &liwork, &info),
      info);

  cudaFreeHost(h_R);
  magma_free_cpu(iwork);

  return EXIT_SUCCESS;
}

template<class T, magma_int_t (*ptr_getrf)(magma_int_t m, magma_int_t N, T *dA,
                                           magma_int_t ldda, magma_int_t *ipiv,
                                           magma_int_t *info),
    magma_int_t (*ptr_getri)(magma_int_t N, T *dA, magma_int_t ldda,
                             magma_int_t *ipiv, T *dwork, magma_int_t lwork,
                             magma_int_t *info),
    magma_int_t (*ptr_get_getri_nb)(magma_int_t N)>
int carma_getri_gen(long N, T *d_iA) {
  magma_int_t *ipiv;
  magma_malloc_cpu((void**) &ipiv, N * sizeof(magma_int_t));

  T *dwork;
  magma_int_t lda = N;
  magma_int_t ldda = ((N + 31) / 32) * 32;
  magma_int_t ldwork = N * ptr_get_getri_nb(N);
  if (lda * ldda != N * N)
    cerr << "Please provide arrays of proper size\n";

  magma_int_t info = 0;
  cudaMalloc((void**) &dwork, ldwork * sizeof(T));
  CHECK_MAGMA(ptr_getrf(N, N, d_iA, ldda, ipiv, &info), info);
  CHECK_MAGMA(ptr_getri(N, d_iA, ldda, ipiv, dwork, ldwork, &info), info);

  cudaFree(dwork);
  magma_free_cpu(ipiv);

  return EXIT_SUCCESS;
}

template<class T, void (*ptr_getrf)(const magma_int_t *m, const magma_int_t *n,
                                    T *h_A, const magma_int_t *lda,
                                    magma_int_t *ipiv, magma_int_t *info),
    void (*ptr_getri)(const magma_int_t *n, T *A, const magma_int_t *lda,
                      const magma_int_t *ipiv, T *work,
                      const magma_int_t *lwork, magma_int_t *info),
    magma_int_t (*ptr_get_getri_nb)(magma_int_t N)>
int carma_getri_cpu_gen(magma_int_t N, T *h_A) {
  magma_int_t *ipiv;
  magma_malloc_cpu((void**) &ipiv, N * sizeof(magma_int_t));

  T *work;
  magma_int_t lda = N;
  magma_int_t lwork = N * ptr_get_getri_nb(N);
  magma_malloc_cpu((void**) &work, lwork * sizeof(T));

  magma_int_t info = 0;
  CHECK_MAGMA(ptr_getrf(&N, &N, h_A, &lda, ipiv, &info), info);
  CHECK_MAGMA(ptr_getri(&N, h_A, &lda, ipiv, work, &lwork, &info), info);

  magma_free_cpu(work);
  magma_free_cpu(ipiv);
  return EXIT_SUCCESS;
}

#if (MAGMA_VERSION_MAJOR == 1) && (MAGMA_VERSION_MINOR <= 5)
template<class T,
magma_int_t (*ptr_potrf)(magma_int_t num_gpus, magma_uplo_t uplo,
    magma_int_t N, T **d_A, magma_int_t ldda, magma_int_t *info),
magma_int_t (*ptr_potri)(magma_uplo_t uplo, magma_int_t N, T *d_A,
    magma_int_t ldda, magma_int_t *info),
magma_int_t (*ptr_get_potrf_nb)(magma_int_t m),
void (*ptr_setmatrix_1D_col_bcyclic)(magma_int_t m, magma_int_t N,
    const T *hA, magma_int_t lda, T *dA[], magma_int_t ldda,
    magma_int_t num_gpus, magma_int_t nb),
void (*ptr_getmatrix_1D_col_bcyclic)(magma_int_t m, magma_int_t N, T *dA[],
    magma_int_t ldda, T *hA, magma_int_t lda, magma_int_t num_gpus,
    magma_int_t nb)>
#else
template<class T, magma_int_t (*ptr_potrf)(magma_int_t num_gpus,
                                           magma_uplo_t uplo, magma_int_t N,
                                           T **d_A, magma_int_t ldda,
                                           magma_int_t *info),
    magma_int_t (*ptr_potri)(magma_uplo_t uplo, magma_int_t N, T *d_A,
                             magma_int_t ldda, magma_int_t *info),
    magma_int_t (*ptr_get_potrf_nb)(magma_int_t m),
    void (*ptr_setmatrix_1D_col_bcyclic)(magma_int_t m, magma_int_t N,
                                         const T *hA, magma_int_t lda, T *dA[],
                                         magma_int_t ldda, magma_int_t num_gpus,
                                         magma_int_t nb),
    void (*ptr_getmatrix_1D_col_bcyclic)(magma_int_t m, magma_int_t N,
                                         const T * const dA[], magma_int_t ldda,
                                         T *hA, magma_int_t lda,
                                         magma_int_t num_gpus, magma_int_t nb)>
#endif
int carma_potri_m_gen(long num_gpus, T *h_A, T *d_iA, long N,
                      carma_device *device) {
  magma_int_t nb = ptr_get_potrf_nb(N);
  magma_int_t ldda = ((N + nb - 1) / nb) * nb;
  magma_int_t lda = N;

  T *d_lA[num_gpus];
  magma_int_t max_size = nb * (1 + N / (nb * num_gpus)) * nb
                         * ((N + nb - 1) / nb);
  for (long dev = 0; dev < num_gpus; dev++) {
    magma_setdevice(dev);
    cudaMalloc(&d_lA[dev], max_size * sizeof(T));
  }
  ptr_setmatrix_1D_col_bcyclic(N, N, h_A, lda, d_lA, ldda, num_gpus, nb);

  magma_int_t info = 0;
  CHECK_MAGMA(ptr_potrf(num_gpus, magma_uplo_const('L'), N, d_lA, N, &info),
              info);
  CHECK_MAGMA(
      ptr_getmatrix_1D_col_bcyclic(N, N, d_lA, ldda, h_A, lda, num_gpus, nb),
      info);

  cudaMemcpy(h_A, d_iA, N * N * sizeof(T), cudaMemcpyHostToDevice);
  //d_iA->host2device(h_A);
  CHECK_MAGMA(ptr_potri(magma_uplo_const('L'), N, d_iA, N, &info), info);

  fill_sym_matrix<T>('L', d_iA, N, N * N, device);
  return EXIT_SUCCESS;
}

template<class T, void (*ptr_potrf)(const char *uplo, const magma_int_t *n,
                                    T *A, const magma_int_t *lda,
                                    magma_int_t *info), void (*ptr_potri)(
    const char *uplo, const magma_int_t *n, T *A, const magma_int_t *lda,
    magma_int_t *info), void (*ptr_copy)(const magma_int_t *n, const T *x,
                                         const magma_int_t *incx, T *y,
                                         const magma_int_t *incy)>
int carma_potri_cpu_gen(magma_int_t N, T *h_A) {
  magma_int_t info = 0;
  char uplo = 'L';
  magma_int_t one = 1;
  CHECK_MAGMA(ptr_potrf(&uplo, &N, h_A, &N, &info), info)
  CHECK_MAGMA(ptr_potri(&uplo, &N, h_A, &N, &info), info)
  // Copy lower part to the upper
  for (magma_int_t row = 1; row < N; row++) {
    ptr_copy(&row, h_A + row, &N, h_A + row * N, &one);
  }
  return EXIT_SUCCESS;
}

#define TEST_USE_MAGMA(...) __VA_ARGS__

#else
#warning "MAGMA will not be used"

#define TEST_USE_MAGMA(...) \
  MAGMA_TRACE("!!!!!! MAGMA not compiled !!!!!!\n"); \
  return EXIT_FAILURE;
#endif

/*
 * API functions
 */

int carma_disabled() {
  TEST_USE_MAGMA(return EXIT_SUCCESS);
}

template<>
int carma_syevd<float>(char jobz, long N, float *mat, float *eigenvals) {
  TEST_USE_MAGMA(return carma_syevd_gen<float,
                 magma_ssyevd_gpu>(magma_vec_const(jobz),N,N,mat,eigenvals));
}
template<>
int carma_syevd<double>(char jobz, long N, double *mat, double *eigenvals) {
  TEST_USE_MAGMA(
      return carma_syevd_gen<double,
      magma_dsyevd_gpu>(magma_vec_const(jobz), N, N, mat, eigenvals));
}

template<>
int carma_syevd<float, 1>(char jobz, caObjS *mat,
                          carma_host_obj<float> *eigenvals) {
  long N = mat->getDims(1);
  if (N != mat->getDims(2)) {
    cerr << "Matrix must be symmetric\n";
    return EXIT_FAILURE;
  }
  TEST_USE_MAGMA(
      return carma_syevd_gen<float,
      magma_ssyevd_gpu>(magma_vec_const(jobz), N, N, *mat, *eigenvals));
}
template<>
int carma_syevd<double, 1>(char jobz, caObjD *mat,
                           carma_host_obj<double> *eigenvals) {
  long N = mat->getDims(1);
  if (N != mat->getDims(2)) {
    cerr << "Matrix must be symmetric\n";
    return EXIT_FAILURE;
  }
  TEST_USE_MAGMA(
      return carma_syevd_gen<double,
      magma_dsyevd_gpu>(magma_vec_const(jobz), N, N, *mat, *eigenvals));
}
/*
 template<> int carma_syevd<float, 2>(char jobz, caObjS *mat,
 carma_host_obj<float> *eigenvals) {
 MAGMA_TRACE("carma_syevd: this method not implemented !\n");
 return EXIT_FAILURE;
 }
 template<> int carma_syevd<double, 2>(char jobz, caObjD *mat,
 carma_host_obj<double> *eigenvals) {
 return carma_syevd<double, magma_dsyevd_gpu_magmablas>(magma_vec_const(jobz), mat, eigenvals);
 }
 template<> int carma_syevd<float, 3>(char jobz, caObjS *mat,
 carma_host_obj<float> *eigenvals) {
 MAGMA_TRACE("carma_syevd: this method not implemented !\n");
 return EXIT_FAILURE;
 }
 template<> int carma_syevd<double, 3>(char jobz, caObjD *mat,
 carma_host_obj<double> *eigenvals) {
 return carma_syevd<double, magma_dsyevd_gpu_kblas>(magma_vec_const(jobz), mat, eigenvals);
 }
 */

template<>
int carma_syevd_m<float>(long ngpu, char jobz, long N, float *mat,
                         float *eigenvals) {
  TEST_USE_MAGMA(
      magma_int_t N_=N; return carma_syevd_m_gen<float,
      magma_ssyevd_m>(ngpu, magma_uplo_const(jobz), N_, mat, eigenvals));
}
template<>
int carma_syevd_m<double>(long ngpu, char jobz, long N, double *mat,
                          double *eigenvals) {
  TEST_USE_MAGMA(
      magma_int_t N_=N; return carma_syevd_m_gen<double,
      magma_dsyevd_m>(ngpu, magma_uplo_const(jobz), N_, mat, eigenvals));
}

template<>
int carma_syevd_m<float>(long ngpu, char jobz, carma_host_obj<float> *mat,
                         carma_host_obj<float> *eigenvals) {
  long N = mat->getDims(1);
  if (N != mat->getDims(2)) {
    cerr << "Matrix must be symmetric\n";
    return EXIT_FAILURE;
  }
  TEST_USE_MAGMA(
      return carma_syevd_m_gen<float,
      magma_ssyevd_m>(ngpu, magma_uplo_const(jobz), N, *mat, *eigenvals));
}
template<>
int carma_syevd_m<double>(long ngpu, char jobz, carma_host_obj<double> *mat,
                          carma_host_obj<double> *eigenvals) {
  long N = mat->getDims(1);
  if (N != mat->getDims(2)) {
    cerr << "Matrix must be symmetric\n";
    return EXIT_FAILURE;
  }
  TEST_USE_MAGMA(
      return carma_syevd_m_gen<double,
      magma_dsyevd_m>(magma_int_t(ngpu), magma_uplo_const(jobz), N, *mat, *eigenvals));
}

template<>
int carma_syevd_m<float>(long ngpu, char jobz, carma_host_obj<float> *mat,
                         carma_host_obj<float> *eigenvals,
                         carma_host_obj<float> *U) {
  TEST_USE_MAGMA(
      /* Initialize the matrix */
      magma_int_t N = mat->getDims(1); if (N != mat->getDims(2)) { cerr << "Matrix must be symmetric\n"; return EXIT_FAILURE; }

      magma_int_t lda = N; lapackf77_slacpy(MagmaUpperLowerStr, &N, &N, *mat, &lda, *U, &lda);
      return carma_syevd_m_gen<float, magma_ssyevd_m>(ngpu, magma_uplo_const(jobz), N, *U, *eigenvals));
}
template<>
int carma_syevd_m<double>(long ngpu, char jobz, carma_host_obj<double> *mat,
                          carma_host_obj<double> *eigenvals,
                          carma_host_obj<double> *U) {
  TEST_USE_MAGMA(
      /* Initialize the matrix */
      magma_int_t N = mat->getDims(1); if (N != mat->getDims(2)) { cerr << "Matrix must be symmetric\n"; return EXIT_FAILURE; }
      magma_int_t lda = N; lapackf77_dlacpy(MagmaUpperLowerStr, &N, &N, *mat, &lda, *U, &lda);
      return carma_syevd_m_gen<double,magma_dsyevd_m>(ngpu, magma_uplo_const(jobz), N, *U, *eigenvals));
}

template<>
int carma_svd<float>(caObjS *imat, caObjS *eigenvals, caObjS *mod2act,
                     caObjS *mes2mod) {
  //TODO: carma_svd
  MAGMA_TRACE("carma_svd not implemented on device object! \n");
  return EXIT_FAILURE;
  //return carma_gesvd<float>(mat, eigenvals, U, magma_sgesvd);
}
template<>
int carma_svd<double>(caObjD *imat, caObjD *eigenvals, caObjD *mod2act,
                      caObjD *mes2mod) {
  //TODO: carma_svd
  MAGMA_TRACE("carma_svd not implemented on device object! \n");
  return EXIT_FAILURE;
  //return carma_gesvd<double>(mat, eigenvals, U, magma_dgesvd);
}

template<>
int carma_svd_cpu<float>(carma_host_obj<float> *mat,
                         carma_host_obj<float> *eigenvals,
                         carma_host_obj<float> *mod2act,
                         carma_host_obj<float> *mes2mod) {
  TEST_USE_MAGMA(return carma_svd_cpu_gen<float,
                 magma_sgesvd>(mat, eigenvals, mod2act, mes2mod));
}
template<>
int carma_svd_cpu<double>(carma_host_obj<double> *mat,
                          carma_host_obj<double> *eigenvals,
                          carma_host_obj<double> *mod2act,
                          carma_host_obj<double> *mes2mod) {
  TEST_USE_MAGMA(return carma_svd_cpu_gen<double,
                 magma_dgesvd>(mat, eigenvals, mod2act, mes2mod));
}

template<>
int carma_potri<float>(carma_obj<float> *d_iA) {
  const long *dims = d_iA->getDims();
  long N = dims[1];
  if (N != dims[2]) {
    MAGMA_TRACE("carma_potri : non square and positive-definite matrix\n");
    return EXIT_FAILURE;
  }
  TEST_USE_MAGMA(
      return carma_potri_gen<float,
      magma_spotrf_gpu,
      magma_spotri_gpu>(N, N, *d_iA, d_iA->getContext()->get_device(d_iA->getDevice())));
}
template<>
int carma_potri<double>(carma_obj<double> *d_iA) {
  const long *dims = d_iA->getDims();
  long N = dims[1];
  if (N != dims[2]) {
    MAGMA_TRACE("carma_potri : non square and positive-definite matrix\n");
    return EXIT_FAILURE;
  }
  TEST_USE_MAGMA(
      return carma_potri_gen<double,
      magma_dpotrf_gpu,
      magma_dpotri_gpu>(N, N, *d_iA, d_iA->getContext()->get_device(d_iA->getDevice())));
}

template<>
int carma_potri_m<float>(long num_gpus, carma_host_obj<float> *h_A,
                         carma_obj<float> *d_iA) {
  long N = h_A->getDims(1);
  if (N != h_A->getDims(2)) {
    cerr << "Matrix must be square and positive-definite\n";
    return EXIT_FAILURE;
  }
  TEST_USE_MAGMA(
      return carma_potri_m_gen<float,
      magma_spotrf_mgpu,
      magma_spotri_gpu,
      magma_get_spotrf_nb,
      magma_ssetmatrix_1D_col_bcyclic,
      magma_sgetmatrix_1D_col_bcyclic>( num_gpus, *h_A, *d_iA, N, d_iA->getContext()->get_device(d_iA->getDevice())));
}
template<>
int carma_potri_m<double>(long num_gpus, carma_host_obj<double> *h_A,
                          carma_obj<double> *d_iA) {
  long N = h_A->getDims(1);
  if (N != h_A->getDims(2)) {
    cerr << "Matrix must be square and positive-definite\n";
    return EXIT_FAILURE;
  }

  TEST_USE_MAGMA(
      return carma_potri_m_gen<double,
      magma_dpotrf_mgpu,
      magma_dpotri_gpu,
      magma_get_dpotrf_nb,
      magma_dsetmatrix_1D_col_bcyclic,
      magma_dgetmatrix_1D_col_bcyclic>( num_gpus, *h_A, *d_iA, N, d_iA->getContext()->get_device(d_iA->getDevice())));
}

template<>
int carma_getri<float>(carma_obj<float> *d_iA) {
  long N = d_iA->getDims(1);
  if (N != d_iA->getDims(2)) {
    cerr << "Matrix must be square\n";
    return EXIT_FAILURE;
  }

  TEST_USE_MAGMA(return carma_getri_gen<float, magma_sgetrf_gpu,
                 magma_sgetri_gpu, magma_get_sgetri_nb>(N, *d_iA));
}
template<>
int carma_getri<double>(carma_obj<double> *d_iA) {
  long N = d_iA->getDims(1);
  if (N != d_iA->getDims(2)) {
    cerr << "Matrix must be square\n";
    return EXIT_FAILURE;
  }
  TEST_USE_MAGMA(return carma_getri_gen<double, magma_dgetrf_gpu,
                 magma_dgetri_gpu, magma_get_dgetri_nb>(N, *d_iA));
}

template<>
int carma_potri_cpu<float>(carma_host_obj<float> *h_A) {
  long N = h_A->getDims(1);
  if (N != h_A->getDims(2)) {
    cerr << "Matrix must be symmetric\n";
    return EXIT_FAILURE;
  }
  TEST_USE_MAGMA(return carma_potri_cpu_gen<float, lapackf77_spotrf,
                 lapackf77_spotri, blasf77_scopy>(N, *h_A));
}
template<>
int carma_potri_cpu<double>(carma_host_obj<double> *h_A) {
  long N = h_A->getDims(1);
  if (N != h_A->getDims(2)) {
    cerr << "Matrix must be symmetric\n";
    return EXIT_FAILURE;
  }
  TEST_USE_MAGMA(return carma_potri_cpu_gen<double, lapackf77_dpotrf,
                 lapackf77_dpotri, blasf77_dcopy>(N, *h_A));
}

template<>
int carma_potri_cpu<float>(long N, float *h_A) {
  TEST_USE_MAGMA(return carma_potri_cpu_gen<float, lapackf77_spotrf,
                 lapackf77_spotri, blasf77_scopy>(N, h_A));
}

template<>
int carma_potri_cpu<double>(long N, double *h_A) {
  TEST_USE_MAGMA(return carma_potri_cpu_gen<double, lapackf77_dpotrf,
                 lapackf77_dpotri, blasf77_dcopy>(N, h_A));
}

template<>
int carma_getri_cpu<float>(carma_host_obj<float> *h_A) {
  long N = h_A->getDims(1);
  if (N != h_A->getDims(2)) {
    cerr << "Matrix must be symmetric\n";
    return EXIT_FAILURE;
  }

  TEST_USE_MAGMA(return carma_getri_cpu_gen<float, lapackf77_sgetrf,
                 lapackf77_sgetri, magma_get_dgetri_nb>(N, *h_A));
}
template<>
int carma_getri_cpu<double>(carma_host_obj<double> *h_A) {
  long N = h_A->getDims(1);
  if (N != h_A->getDims(2)) {
    cerr << "Matrix must be symmetric\n";
    return EXIT_FAILURE;
  }
  TEST_USE_MAGMA(return carma_getri_cpu_gen<double, lapackf77_dgetrf,
                 lapackf77_dgetri, magma_get_sgetri_nb>(N, *h_A));
}

template<>
int carma_getri_cpu<float>(long N, float *h_A) {
  TEST_USE_MAGMA(return carma_getri_cpu_gen<float, lapackf77_sgetrf,
                 lapackf77_sgetri, magma_get_dgetri_nb>(N, h_A));
}

template<>
int carma_getri_cpu<double>(long N, double *h_A) {
  TEST_USE_MAGMA(return carma_getri_cpu_gen<double, lapackf77_dgetrf,
                 lapackf77_dgetri, magma_get_sgetri_nb>(N, h_A));
}

template<>
int carma_syevd_cpu<float>(char jobz, carma_host_obj<float> *h_A,
                           carma_host_obj<float> *eigenvals) {
  long N = h_A->getDims(1);
  if (N != h_A->getDims(2)) {
    cerr << "Matrix must be symmetric\n";
    return EXIT_FAILURE;
  }
  TEST_USE_MAGMA(
      return carma_syevd_cpu_gen<float,
      lapackf77_ssyevd>(magma_vec_const(jobz), N, N, *h_A, *eigenvals));
}
template<>
int carma_syevd_cpu<double>(char jobz, carma_host_obj<double> *h_A,
                            carma_host_obj<double> *eigenvals) {
  long N = h_A->getDims(1);
  if (N != h_A->getDims(2)) {
    cerr << "Matrix must be symmetric\n";
    return EXIT_FAILURE;
  }
  TEST_USE_MAGMA(
      return carma_syevd_cpu_gen<double,
      lapackf77_dsyevd>(magma_vec_const(jobz), N, N, *h_A, *eigenvals));
}

template<>
int carma_syevd_cpu<float>(char jobz, long N, float *h_A, float *eigenvals) {
  TEST_USE_MAGMA(
      magma_int_t N_ = N; return carma_syevd_cpu_gen<float,
      lapackf77_ssyevd>(magma_vec_const(jobz), N_, N_, h_A, eigenvals));
}
template<>
int carma_syevd_cpu<double>(char jobz, long N, double *h_A, double *eigenvals) {
  TEST_USE_MAGMA(
      magma_int_t N_ = N; return carma_syevd_cpu_gen<double,
      lapackf77_dsyevd>(magma_vec_const(jobz), N_, N_, h_A, eigenvals));
}

template<>
int carma_axpy_cpu<float>(long N, float alpha, float *h_X, long incX,
                          float *h_Y, long incY) {
  TEST_USE_MAGMA(
      magma_int_t tmp_N = N; magma_int_t tmp_incX = incX; magma_int_t tmp_incY = incY;
      blasf77_saxpy(&tmp_N, &alpha, h_X, &tmp_incX, h_Y, &tmp_incY));
  return EXIT_SUCCESS;
}
template<>
int carma_axpy_cpu<double>(long N, double alpha, double *h_X, long incX,
                           double *h_Y, long incY) {
  TEST_USE_MAGMA(
      magma_int_t tmp_N = N; magma_int_t tmp_incX = incX; magma_int_t tmp_incY = incY;
      blasf77_daxpy(&tmp_N, &alpha, h_X, &tmp_incX, h_Y, &tmp_incY));

  return EXIT_SUCCESS;
}

template<>
int carma_gemm_cpu<float>(char transa, char transb, long m, long n, long k,
                          float alpha, float *A, long lda, float *B, long ldb,
                          float beta, float *C, long ldc) {
  TEST_USE_MAGMA(
      magma_int_t tmp_m = m; magma_int_t tmp_n = n; magma_int_t tmp_k = k;
      magma_int_t tmp_lda = lda; magma_int_t tmp_ldb = ldb; magma_int_t tmp_ldc = ldc;
      blasf77_sgemm(&transa, &transb, &tmp_m, &tmp_n, &tmp_k, &alpha, A, &tmp_lda, B, &tmp_ldb, &beta, C, &tmp_ldc));
  return EXIT_SUCCESS;
}
template<>
int carma_gemm_cpu<double>(char transa, char transb, long m, long n, long k,
                           double alpha, double *A, long lda, double *B,
                           long ldb, double beta, double *C, long ldc) {
  TEST_USE_MAGMA(
      magma_int_t tmp_m = m; magma_int_t tmp_n = n; magma_int_t tmp_k = k;
      magma_int_t tmp_lda = lda; magma_int_t tmp_ldb = ldb; magma_int_t tmp_ldc = ldc;
      blasf77_dgemm(&transa, &transb, &tmp_m, &tmp_n, &tmp_k, &alpha, A, &tmp_lda, B, &tmp_ldb, &beta, C, &tmp_ldc));
  return EXIT_SUCCESS;
}

