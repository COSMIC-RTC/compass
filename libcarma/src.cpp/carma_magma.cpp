#include <carma_obj.h>
#include <carma_host_obj.h>

#define MAGMA_TRACE(fmt, args...) fprintf(stderr, "%s:%d Warning: " fmt, __FILE__, __LINE__, ## args)

#ifdef USE_MAGMA
#include "magma.h"
#include "magma_lapack.h"

#define CHECK_MAGMA(fct, info) fct; \
							   if (info != 0){ \
								   printf("%s@%d magma returned error %d: %s.\n", __FILE__, __LINE__, \
										   (int) info, magma_strerror(info)); \
								   return EXIT_FAILURE; \
							   }

#ifndef max
#define max(a,b)  (((a)<(b))?(b):(a))
#endif
#ifndef min
#define min(a,b)  (((a)>(b))?(b):(a))
#endif

template<class T>
  T*
  create_padded_data(magma_int_t N, const magma_int_t size, magma_int_t &ldda,
      T* d_src_data) {
    ldda = ((N + (size - 1)) / size) * size;
//	MAGMA_TRACE(
//			"the dimension of the matrix (N=%ld) is not a multiple of %ld -> padding leading dimension to %ld\n",
//			(long)N, (long)size, (long)ldda);
    T *d_padded_data;
    magma_malloc((void**) &d_padded_data, N * ldda * sizeof(T));
    cudaMemcpy2D(d_padded_data, ldda * sizeof(T), d_src_data, N * sizeof(T),
        N * sizeof(T), N, cudaMemcpyDeviceToDevice);
    return d_padded_data;
  }
template<class T>
  void
  destroy_padded_data(magma_int_t N, magma_int_t ldda, T* d_padded_data,
      T* d_src_data) {
    cudaMemcpy2D(d_src_data, N * sizeof(T), d_padded_data, ldda * sizeof(T),
        N * sizeof(T), N, cudaMemcpyDeviceToDevice);
    magma_free(d_padded_data);
  }

template<class T>
  int
  carma_syevd_gen(magma_vec_t jobz, magma_int_t N, magma_int_t ldda, T *d_mat,
      T *h_eigenvals,
      magma_int_t
      (*ptr_syevd_gpu)(magma_vec_t jobz, magma_uplo_t uplo, magma_int_t n,
          T *da, magma_int_t ldda, T *w, T *wa, magma_int_t ldwa, T *work,
          magma_int_t lwork, magma_int_t *iwork, magma_int_t liwork,
          magma_int_t *info)) {
    magma_int_t *iwork;
//	if (ldda & 31) {
//		cerr << "Leading dimension of the matrix must be a multiple of 32\n";
//		return EXIT_FAILURE;
//	}

    magma_int_t info = 0, lwork, liwork, lda = N, aux_iwork[1];
    T aux_work[1];
    T *h_R, *h_work;

    /* Query for workspace sizes */
    CHECK_MAGMA(
        ptr_syevd_gpu(jobz, 'L', N, NULL, ldda, NULL, NULL, lda, aux_work, -1, aux_iwork, -1, &info),
        info);

    lwork = (magma_int_t) aux_work[0];
    liwork = aux_iwork[0];
    magma_malloc_cpu((void**) &iwork, (liwork) * sizeof(magma_int_t));

    cudaMallocHost((void**) &h_R, (N * lda) * sizeof(T));
    cudaMallocHost((void**) &h_work, (lwork) * sizeof(T));

    CHECK_MAGMA(
        ptr_syevd_gpu(jobz, 'L', N, d_mat, ldda, h_eigenvals, h_R, lda, h_work,
            lwork, iwork, liwork, &info), info);

    cudaFreeHost(h_R);
    cudaFreeHost(h_work);
    free(iwork);

    return EXIT_SUCCESS;
  }

template<class T>
  int
  carma_syevd(char jobz, carma_obj<T> *mat, carma_host_obj<T> *eigenvals,
      magma_int_t
      (*ptr_syevd_gpu)(magma_vec_t jobz, magma_uplo_t uplo, magma_int_t n,
          T *da, magma_int_t ldda, T *w, T *wa, magma_int_t ldwa, T *work,
          magma_int_t lwork, magma_int_t *iwork, magma_int_t liwork,
          magma_int_t *info)) {
    cerr << "carma_syevd: this method not implemented !\n";
    return EXIT_FAILURE;
  }
template<>
  int
  carma_syevd<float>(char jobz, caObjS *mat, carma_host_obj<float> *eigenvals,
      magma_int_t
      (*ptr_syevd_gpu)(magma_vec_t jobz, magma_uplo_t uplo, magma_int_t n,
          float *da, magma_int_t ldda, float *w, float *wa, magma_int_t ldwa,
          float *work, magma_int_t lwork, magma_int_t *iwork,
          magma_int_t liwork, magma_int_t *info)) {
    magma_int_t N = mat->getDims(1);
    if (N != mat->getDims(2)) {
      cerr << "Matrix must be symmetric\n";
      return EXIT_FAILURE;
    }
    /*
     if (N & 31) { //if dimension of mat is not a multiple of 32
     magma_int_t ldda = 0;
     float *d_padded_data = create_padded_data<float>(N, 32, ldda, *mat);
     int err = carma_syevd_gen<float>(jobz, N, ldda, d_padded_data,
     eigenvals, ptr_syevd_gpu);
     destroy_padded_data<float>(N, ldda, d_padded_data, *mat);
     return err;
     }
     */
    return carma_syevd_gen<float>(jobz, N, N, *mat, *eigenvals, ptr_syevd_gpu);
  }
template<>
  int
  carma_syevd<double>(char jobz, caObjD *mat, carma_host_obj<double> *eigenvals,
      magma_int_t
      (*ptr_syevd_gpu)(magma_vec_t jobz, magma_uplo_t uplo, magma_int_t n,
          double *da, magma_int_t ldda, double *w, double *wa, magma_int_t ldwa,
          double *work, magma_int_t lwork, magma_int_t *iwork,
          magma_int_t liwork, magma_int_t *info)) {
    magma_int_t N = mat->getDims(1);
    if (N != mat->getDims(2)) {
      cerr << "Matrix must be symmetric\n";
      return EXIT_FAILURE;
    }
    /*
     if (N & 31) { //if dimension of mat is not a multiple of 32
     magma_int_t ldda = 0;
     double *d_padded_data = create_padded_data<double>(N, 32, ldda, *mat);
     int err = carma_syevd_gen<double>(jobz, N, ldda, d_padded_data,
     eigenvals, ptr_syevd_gpu);
     destroy_padded_data<double>(N, ldda, d_padded_data, *mat);
     return err;
     }
     */
    return carma_syevd_gen<double>(jobz, N, N, *mat, *eigenvals, ptr_syevd_gpu);
  }

template<class T>
  int
  carma_syevd(char jobz, magma_int_t N, T *mat, T *eigenvals) {
    cerr << "carma_syevd: this method not implemented !\n";
    return EXIT_FAILURE;
  }
template<>
  int
  carma_syevd<float>(char jobz, magma_int_t N, float *mat, float *eigenvals) {
    return carma_syevd_gen<float>(jobz, N, N, mat, eigenvals, magma_ssyevd_gpu);
  }
template<>
  int
  carma_syevd<double>(char jobz, magma_int_t N, double *mat,
      double *eigenvals) {
    return carma_syevd_gen<double>(jobz, N, N, mat, eigenvals, magma_dsyevd_gpu);
  }

template<class T, int method>
  int
  carma_syevd(char jobz, carma_obj<T> *mat, carma_host_obj<T> *eigenvals) {
    MAGMA_TRACE("carma_syevd: this method not implemented !\n");
    return EXIT_FAILURE;
  }
template<>
  int
  carma_syevd<float, 1>(char jobz, caObjS *mat,
      carma_host_obj<float> *eigenvals) {
    return carma_syevd_gen<float>(jobz, N, N, mat, eigenvals, magma_ssyevd_gpu);
  }
template<>
  int
  carma_syevd<double, 1>(char jobz, caObjD *mat,
      carma_host_obj<double> *eigenvals) {
    return carma_syevd_gen<double>(jobz, N, N, mat, eigenvals, magma_dsyevd_gpu);
  }
/*
 template<> int carma_syevd<float, 2>(char jobz, caObjS *mat,
 carma_host_obj<float> *eigenvals) {
 MAGMA_TRACE("carma_syevd: this method not implemented !\n");
 return EXIT_FAILURE;
 }
 template<> int carma_syevd<double, 2>(char jobz, caObjD *mat,
 carma_host_obj<double> *eigenvals) {
 return carma_syevd<double>(jobz, mat, eigenvals, magma_dsyevd_gpu_magmablas);
 }
 template<> int carma_syevd<float, 3>(char jobz, caObjS *mat,
 carma_host_obj<float> *eigenvals) {
 MAGMA_TRACE("carma_syevd: this method not implemented !\n");
 return EXIT_FAILURE;
 }
 template<> int carma_syevd<double, 3>(char jobz, caObjD *mat,
 carma_host_obj<double> *eigenvals) {
 return carma_syevd<double>(jobz, mat, eigenvals, magma_dsyevd_gpu_kblas);
 }
 */
template<class T>
  int
  carma_syevd_m_gen(magma_int_t ngpu, magma_uplo_t jobz, magma_int_t N, T *mat,
      T *eigenvals,
      magma_int_t
      (*ptr_syevd_m)(magma_int_t nrgpu, magma_vec_t jobz, magma_uplo_t uplo,
          magma_int_t n, T *a, magma_int_t lda, T *w, T *work,
          magma_int_t lwork, magma_int_t *iwork, magma_int_t liwork,
          magma_int_t *info)) {
    magma_int_t *iwork;
    magma_int_t info = 0, lwork, liwork, lda = N, aux_iwork[1];
    T aux_work[1];
    T *h_work;

    /* Query for workspace sizes */
    CHECK_MAGMA(
        ptr_syevd_m(ngpu, jobz, 'L', N, NULL, lda, NULL, aux_work, -1, aux_iwork, -1, &info),
        info);

    //fprintf(stderr, "%s@%d : I'm here\n", __FILE__, __LINE__);
    lwork = (magma_int_t) aux_work[0];
    liwork = aux_iwork[0];

    magma_malloc_cpu((void**) &iwork, (liwork) * sizeof(magma_int_t));
    //cudaMallocHost((void**) &h_R, (N * lda) * sizeof(T));
    cudaMallocHost((void**) &h_work, (lwork) * sizeof(T));

    CHECK_MAGMA(
        ptr_syevd_m(ngpu, jobz, 'L', N, mat, lda, eigenvals, h_work, lwork,
            iwork, liwork, &info), info);

    //cudaFreeHost(h_R);
    cudaFreeHost(h_work);
    free(iwork);

    return EXIT_SUCCESS;
  }

template<class T>
  int
  carma_syevd_m(long ngpu, char jobz, magma_int_t N, T *mat, T *eigenvals) {
    MAGMA_TRACE("carma_syevd_m not implemented this this type! \n");
    return EXIT_FAILURE;
  }
template<>
  int
  carma_syevd_m<float>(long ngpu, char jobz, magma_int_t N, float *mat,
      float *eigenvals) {
    return carma_syevd_m_gen<float>(ngpu, jobz, N, mat, eigenvals,
        magma_ssyevd_m);
  }
template<>
  int
  carma_syevd_m<double>(long ngpu, char jobz, magma_int_t N, double *mat,
      double *eigenvals) {
    return carma_syevd_m_gen<double>(ngpu, jobz, N, mat, eigenvals,
        magma_dsyevd_m);
  }

template<class T>
  int
  carma_syevd_m(long ngpu, char jobz, carma_host_obj<T> *mat,
      carma_host_obj<T> *eigenvals) {
    MAGMA_TRACE("carma_syevd_m not implemented this this type! \n");
    return EXIT_FAILURE;
  }
template<>
  int
  carma_syevd_m<float>(long ngpu, char jobz, carma_host_obj<float> *mat,
      carma_host_obj<float> *eigenvals) {
    magma_int_t N = mat->getDims(1);
    if (N != mat->getDims(2)) {
      cerr << "Matrix must be symmetric\n";
      return EXIT_FAILURE;
    }
    return carma_syevd_m_gen<float>(ngpu, jobz, N, *mat, *eigenvals,
        magma_ssyevd_m);
  }
template<>
  int
  carma_syevd_m<double>(long ngpu, char jobz, carma_host_obj<double> *mat,
      carma_host_obj<double> *eigenvals) {
    magma_int_t N = mat->getDims(1);
    if (N != mat->getDims(2)) {
      cerr << "Matrix must be symmetric\n";
      return EXIT_FAILURE;
    }
    return carma_syevd_m_gen<double>(ngpu, jobz, N, *mat, *eigenvals,
        magma_dsyevd_m);
  }

template<class T>
  int
  carma_syevd_m(long ngpu, char jobz, carma_host_obj<T> *mat,
      carma_host_obj<T> *eigenvals, carma_host_obj<T> *U) {
    MAGMA_TRACE("carma_syevd_m not implemented this this type! \n");
    return EXIT_FAILURE;
  }
template<>
  int
  carma_syevd_m<float>(long ngpu, char jobz, carma_host_obj<float> *mat,
      carma_host_obj<float> *eigenvals, carma_host_obj<float> *U) {
    /* Initialize the matrix */
    magma_int_t N = mat->getDims(1);
    if (N != mat->getDims(2)) {
      cerr << "Matrix must be symmetric\n";
      return EXIT_FAILURE;
    }

    magma_int_t lda = N;
    lapackf77_slacpy(MagmaUpperLowerStr, &N, &N, *mat, &lda, *U, &lda);

    return carma_syevd_m_gen<float>(ngpu, jobz, N, *U, *eigenvals,
        magma_ssyevd_m);
  }
template<>
  int
  carma_syevd_m<double>(long ngpu, char jobz, carma_host_obj<double> *mat,
      carma_host_obj<double> *eigenvals, carma_host_obj<double> *U) {
    /* Initialize the matrix */
    magma_int_t N = mat->getDims(1);
    if (N != mat->getDims(2)) {
      cerr << "Matrix must be symmetric\n";
      return EXIT_FAILURE;
    }
    magma_int_t lda = N;
    lapackf77_dlacpy(MagmaUpperLowerStr, &N, &N, *mat, &lda, *U, &lda);

    return carma_syevd_m_gen<double>(ngpu, jobz, N, *U, *eigenvals,
        magma_dsyevd_m);
  }

template<class T>
  int
  carma_svd(carma_obj<T> *mat, carma_obj<T> *eigenvals, carma_obj<T> *mod2act,
      carma_obj<T> *mes2mod) {
    //TODO: carma_svd
    MAGMA_TRACE("carma_svd not implemented this this type! \n");
    return EXIT_FAILURE;
  }

template<>
  int
  carma_svd<float>(caObjS *imat, caObjS *eigenvals, caObjS *mod2act,
      caObjS *mes2mod) {
    MAGMA_TRACE("carma_svd not implemented on device object! \n");
    return EXIT_FAILURE;
    //return carma_gesvd<float>(mat, eigenvals, U, magma_sgesvd);
  }
template<>
  int
  carma_svd<double>(caObjD *imat, caObjD *eigenvals, caObjD *mod2act,
      caObjD *mes2mod) {
    MAGMA_TRACE("carma_svd not implemented on device object! \n");
    return EXIT_FAILURE;
    //return carma_gesvd<double>(mat, eigenvals, U, magma_dgesvd);
  }

template<class T>
  int
  carma_svd_gen(carma_host_obj<T> *mat, carma_host_obj<T> *eigenvals,
      carma_host_obj<T> *mod2act, carma_host_obj<T> *mes2mod,
      magma_int_t
      (*ptr_gesvd)(magma_vec_t jobu, magma_vec_t jobvt, magma_int_t m,
          magma_int_t n, T *A, magma_int_t lda, T *s, T *U, magma_int_t ldu,
          T *VT, magma_int_t ldvt, T *work, magma_int_t lwork,
          magma_int_t *info)) {
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
        ptr_gesvd('A', 'A', m, n, *tmp, m, *eigenvals, *mes2mod, m, *mod2act, n,
            *h_work, lwork, &info), info);

    delete h_work;
    delete tmp;

    return EXIT_SUCCESS;
  }

template<class T>
  int
  carma_svd(carma_host_obj<T> *mat, carma_host_obj<T> *eigenvals,
      carma_host_obj<T> *mod2act, carma_host_obj<T> *mes2mod) {
    //TODO: carma_svd
    MAGMA_TRACE("carma_svd not implemented with this type! \n");
    return EXIT_FAILURE;
  }
template<>
  int
  carma_svd<float>(carma_host_obj<float> *mat, carma_host_obj<float> *eigenvals,
      carma_host_obj<float> *mod2act, carma_host_obj<float> *mes2mod) {
    return carma_svd_gen<float>(mat, eigenvals, mod2act, mes2mod, magma_sgesvd);
  }
template<>
  int
  carma_svd<double>(carma_host_obj<double> *mat,
      carma_host_obj<double> *eigenvals, carma_host_obj<double> *mod2act,
      carma_host_obj<double> *mes2mod) {
    return carma_svd_gen<double>(mat, eigenvals, mod2act, mes2mod, magma_dgesvd);
  }

template<class T>
  int
  carma_potri_gen(magma_int_t N, magma_int_t ldda, T *d_iA,
      magma_int_t
      (*ptr_potrf)(magma_uplo_t uplo, magma_int_t n, T *d_A, magma_int_t ldda,
          magma_int_t *info),
      magma_int_t
      (*ptr_potri)(magma_uplo_t uplo, magma_int_t n, T *d_A, magma_int_t ldda,
          magma_int_t *info)) {

    magma_int_t info = 0;
    CHECK_MAGMA(ptr_potrf('L', N, d_iA, ldda, &info), info);
    CHECK_MAGMA(ptr_potri('L', N, d_iA, ldda, &info), info);

    fill_sym_matrix<T>('L', d_iA, N, N * ldda);
    return EXIT_SUCCESS;
  }

template<class T>
  int
  carma_potri(carma_obj<T> *d_iA) {
    MAGMA_TRACE("carma_potri not implemented with this type! \n");
    return EXIT_FAILURE;
  }
template<>
  int
  carma_potri<float>(carma_obj<float> *d_iA) {
    const long *dims = d_iA->getDims();
    long N = dims[1];
    if (N != dims[2]) {
      MAGMA_TRACE("carma_potri : non symmetric matrix\n");
      return EXIT_FAILURE;
    }
    return carma_potri_gen<float>(N, N, *d_iA, magma_spotrf_gpu,
        magma_spotri_gpu);
  }
template<>
  int
  carma_potri<double>(carma_obj<double> *d_iA) {
    const long *dims = d_iA->getDims();
    long N = dims[1];
    if (N != dims[2]) {
      MAGMA_TRACE("carma_potri : non symmetric matrix\n");
      return EXIT_FAILURE;
    }
    return carma_potri_gen<double>(N, N, *d_iA, magma_dpotrf_gpu,
        magma_dpotri_gpu);
  }

template<class T>
  int
  carma_potri_m_gen(long num_gpus, T *h_A, T *d_iA, long N,
      magma_int_t
      (*ptr_potrf)(magma_int_t num_gpus, magma_uplo_t uplo, magma_int_t N,
          T **d_A, magma_int_t ldda, magma_int_t *info),
      magma_int_t
      (*ptr_potri)(magma_uplo_t uplo, magma_int_t N, T *d_A, magma_int_t ldda,
          magma_int_t *info), magma_int_t
      (*ptr_get_potrf_nb)(magma_int_t m),
      void
      (*ptr_setmatrix_1D_col_bcyclic)(magma_int_t m, magma_int_t N, const T *hA,
          magma_int_t lda, T *dA[], magma_int_t ldda, magma_int_t num_gpus,
          magma_int_t nb),
      void
      (*ptr_getmatrix_1D_col_bcyclic)(magma_int_t m, magma_int_t N, T *dA[],
          magma_int_t ldda, T *hA, magma_int_t lda, magma_int_t num_gpus,
          magma_int_t nb)) {
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
    CHECK_MAGMA(ptr_potrf(num_gpus, 'L', N, d_lA, N, &info), info);
    CHECK_MAGMA(
        ptr_getmatrix_1D_col_bcyclic(N, N, d_lA, ldda, h_A, lda, num_gpus, nb),
        info);

    cudaMemcpy(h_A, d_iA, N * N * sizeof(T), cudaMemcpyHostToDevice);
    //d_iA->host2device(h_A);
    CHECK_MAGMA(ptr_potri('L', N, d_iA, N, &info), info);

    fill_sym_matrix<T>('L', d_iA, N, N * N);
    return EXIT_SUCCESS;
  }

template<class T>
  int
  carma_potri_m(long num_gpus, carma_host_obj<T> *h_A, carma_obj<T> *d_iA) {
    MAGMA_TRACE("carma_potri_m not implemented with this type!\n");
    return EXIT_FAILURE;
  }

template<>
  int
  carma_potri_m<float>(long num_gpus, carma_host_obj<float> *h_A,
      carma_obj<float> *d_iA) {
    magma_int_t N = h_A->getDims(1);
    if (N != h_A->getDims(2)) {
      cerr << "Matrix must be symmetric\n";
      return EXIT_FAILURE;
    }
    return carma_potri_m_gen<float>(num_gpus, *h_A, *d_iA, N, magma_spotrf_mgpu,
        magma_spotri_gpu, magma_get_spotrf_nb, magma_ssetmatrix_1D_col_bcyclic,
        magma_sgetmatrix_1D_col_bcyclic);
  }
template<>
  int
  carma_potri_m<double>(long num_gpus, carma_host_obj<double> *h_A,
      carma_obj<double> *d_iA) {
    magma_int_t N = h_A->getDims(1);
    if (N != h_A->getDims(2)) {
      cerr << "Matrix must be symmetric\n";
      return EXIT_FAILURE;
    }

    return carma_potri_m_gen<double>(num_gpus, *h_A, *d_iA, N,
        magma_dpotrf_mgpu, magma_dpotri_gpu, magma_get_dpotrf_nb,
        magma_dsetmatrix_1D_col_bcyclic, magma_dgetmatrix_1D_col_bcyclic);
  }

template<class T>
  int
  carma_getri_gen(long N, T *d_iA,
      magma_int_t
      (*ptr_getrf)(magma_int_t m, magma_int_t N, T *dA, magma_int_t ldda,
          magma_int_t *ipiv, magma_int_t *info),
      magma_int_t
      (*ptr_getri)(magma_int_t N, T *dA, magma_int_t ldda, magma_int_t *ipiv,
          T *dwork, magma_int_t lwork, magma_int_t *info), magma_int_t
      (*ptr_get_getri_nb)(magma_int_t N)) {
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
    free(ipiv);

    return EXIT_SUCCESS;
  }

template<class T>
  int
  carma_getri(carma_obj<T> *d_iA) {
    MAGMA_TRACE("carma_getri not implemented with this type! \n");
    return EXIT_FAILURE;
  }
template<>
  int
  carma_getri<float>(carma_obj<float> *d_iA) {
    magma_int_t N = d_iA->getDims(1);
    if (N != d_iA->getDims(2)) {
      cerr << "Matrix must be symmetric\n";
      return EXIT_FAILURE;
    }

    return carma_getri_gen<float>(N, *d_iA, magma_sgetrf_gpu, magma_sgetri_gpu,
        magma_get_sgetri_nb);
  }
template<>
  int
  carma_getri<double>(carma_obj<double> *d_iA) {
    magma_int_t N = d_iA->getDims(1);
    if (N != d_iA->getDims(2)) {
      cerr << "Matrix must be symmetric\n";
      return EXIT_FAILURE;
    }

    return carma_getri_gen<double>(N, *d_iA, magma_dgetrf_gpu, magma_dgetri_gpu,
        magma_get_dgetri_nb);
  }

template<class T>
  int
  carma_potri_cpu_gen(magma_int_t N, T *h_A,
      void
      (*ptr_potrf)(const char *uplo, const magma_int_t *n, T *A,
          const magma_int_t *lda, magma_int_t *info),
      void
      (*ptr_potri)(const char *uplo, const magma_int_t *n, T *A,
          const magma_int_t *lda, magma_int_t *info),
      void
      (*ptr_copy)(const magma_int_t *n, const T *x, const magma_int_t *incx,
          T *y, const magma_int_t *incy)) {
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

template<class T>
  int
  carma_potri_cpu(carma_host_obj<T> *h_A) {
    MAGMA_TRACE("carma_potri not implemented with this type! \n");
    return EXIT_FAILURE;
  }

template<>
  int
  carma_potri_cpu<float>(carma_host_obj<float> *h_A) {
    magma_int_t N = h_A->getDims(1);
    if (N != h_A->getDims(2)) {
      cerr << "Matrix must be symmetric\n";
      return EXIT_FAILURE;
    }

    return carma_potri_cpu_gen<float>(N, *h_A, lapackf77_spotrf,
        lapackf77_spotri, blasf77_scopy);
  }
template<>
  int
  carma_potri_cpu<double>(carma_host_obj<double> *h_A) {
    magma_int_t N = h_A->getDims(1);
    if (N != h_A->getDims(2)) {
      cerr << "Matrix must be symmetric\n";
      return EXIT_FAILURE;
    }
    return carma_potri_cpu_gen<double>(N, *h_A, lapackf77_dpotrf,
        lapackf77_dpotri, blasf77_dcopy);
  }

template<class T>
  int
  carma_potri_cpu(long N, T *h_A) {
    MAGMA_TRACE("carma_potri not implemented with this type! \n");
    return EXIT_FAILURE;
  }

template<>
  int
  carma_potri_cpu<float>(long N, float *h_A) {
    return carma_potri_cpu_gen<float>(N, h_A, lapackf77_spotrf,
        lapackf77_spotri, blasf77_scopy);
  }

template<>
  int
  carma_potri_cpu<double>(long N, double *h_A) {
    return carma_potri_cpu_gen<double>(N, h_A, lapackf77_dpotrf,
        lapackf77_dpotri, blasf77_dcopy);
  }

template<class T>
  int
  carma_getri_cpu_gen(magma_int_t N, T *h_A,
      void
      (*ptr_getrf)(const magma_int_t *m, const magma_int_t *n, T *h_A,
          const magma_int_t *lda, magma_int_t *ipiv, magma_int_t *info),
      void
      (*ptr_getri)(const magma_int_t *n, T *A, const magma_int_t *lda,
          const magma_int_t *ipiv, T *work, const magma_int_t *lwork,
          magma_int_t *info), magma_int_t
      (*ptr_get_getri_nb)(magma_int_t N)) {
    magma_int_t *ipiv;
    magma_malloc_cpu((void**) &ipiv, N * sizeof(magma_int_t));

    T *work;
    magma_int_t lda = N;
    magma_int_t lwork = N * ptr_get_getri_nb(N);
    magma_malloc_cpu((void**) &work, lwork * sizeof(T));

    magma_int_t info = 0;
    CHECK_MAGMA(ptr_getrf(&N, &N, h_A, &lda, ipiv, &info), info);
    CHECK_MAGMA(ptr_getri(&N, h_A, &lda, ipiv, work, &lwork, &info), info);

    free(work);
    free(ipiv);
    return EXIT_SUCCESS;
  }

template<class T>
  int
  carma_getri_cpu(carma_host_obj<T> *h_A) {
    MAGMA_TRACE("carma_potri not implemented with this type! \n");
    return EXIT_FAILURE;
  }

template<>
  int
  carma_getri_cpu<float>(carma_host_obj<float> *h_A) {
    magma_int_t N = h_A->getDims(1);
    if (N != h_A->getDims(2)) {
      cerr << "Matrix must be symmetric\n";
      return EXIT_FAILURE;
    }

    return carma_getri_cpu_gen<float>(N, *h_A, lapackf77_sgetrf,
        lapackf77_sgetri, magma_get_dgetri_nb);
  }
template<>
  int
  carma_getri_cpu<double>(carma_host_obj<double> *h_A) {
    magma_int_t N = h_A->getDims(1);
    if (N != h_A->getDims(2)) {
      cerr << "Matrix must be symmetric\n";
      return EXIT_FAILURE;
    }
    return carma_getri_cpu_gen<double>(N, *h_A, lapackf77_dgetrf,
        lapackf77_dgetri, magma_get_sgetri_nb);
  }

template<class T>
  int
  carma_getri_cpu(long N, T *h_A) {
    MAGMA_TRACE("carma_potri not implemented with this type! \n");
    return EXIT_FAILURE;
  }

template<>
  int
  carma_getri_cpu<float>(long N, float *h_A) {
    return carma_getri_cpu_gen<float>(N, h_A, lapackf77_sgetrf,
        lapackf77_sgetri, magma_get_dgetri_nb);
  }

template<>
  int
  carma_getri_cpu<double>(long N, double *h_A) {
    return carma_getri_cpu_gen<double>(N, h_A, lapackf77_dgetrf,
        lapackf77_dgetri, magma_get_sgetri_nb);
  }

template<class T>
  int
  carma_syevd_cpu_gen(magma_vec_t jobz, magma_int_t N, magma_int_t lda,
      T *d_mat, T *h_eigenvals,
      void
      (*ptr_syevd_cpu)(const char *jobz, const char *uplo, const magma_int_t *n,
          T *da, const magma_int_t *lda, T *w, T *wa, const magma_int_t *ldwa,
          magma_int_t *iwork, const magma_int_t *liwork, magma_int_t *info)) {
    magma_int_t *iwork;

    magma_int_t info = 0, liwork = -1, aux_iwork[1];
    T *h_R;
    DEBUG_TRACE("ICI\n");

    /* Query for workspace sizes */
    char uplo = 'L';
    CHECK_MAGMA(
        ptr_syevd_cpu(&jobz, &uplo, &N, NULL, &lda, NULL, NULL, &lda, aux_iwork, &liwork, &info),
        info);
    DEBUG_TRACE("ICI\n");

    liwork = aux_iwork[0];
    magma_malloc_cpu((void**) &iwork, (liwork) * sizeof(magma_int_t));

    cudaMallocHost((void**) &h_R, (N * lda) * sizeof(T));
    DEBUG_TRACE("ICI\n");

    CHECK_MAGMA(
        ptr_syevd_cpu(&jobz, &uplo, &N, d_mat, &lda, h_eigenvals, h_R, &lda,
            iwork, &liwork, &info), info);

    DEBUG_TRACE("ICI\n");
    cudaFreeHost(h_R);
    free(iwork);
    DEBUG_TRACE("ICI\n");

    return EXIT_SUCCESS;
  }

template<class T>
  int
  carma_syevd_cpu(char jobz, carma_host_obj<T> *h_A,
      carma_host_obj<T> *eigenvals) {
    MAGMA_TRACE("carma_getri not implemented with this type! \n");
    return EXIT_FAILURE;
  }
template<>
  int
  carma_syevd_cpu<float>(char jobz, carma_host_obj<float> *h_A,
      carma_host_obj<float> *eigenvals) {
    magma_int_t N = h_A->getDims(1);
    if (N != h_A->getDims(2)) {
      cerr << "Matrix must be symmetric\n";
      return EXIT_FAILURE;
    }
    return carma_syevd_cpu_gen<float>(jobz, N, N, *h_A, *eigenvals,
        lapackf77_ssyevd);
  }
template<>
  int
  carma_syevd_cpu<double>(char jobz, carma_host_obj<double> *h_A,
      carma_host_obj<double> *eigenvals) {
    magma_int_t N = h_A->getDims(1);
    if (N != h_A->getDims(2)) {
      cerr << "Matrix must be symmetric\n";
      return EXIT_FAILURE;
    }
    return carma_syevd_cpu_gen<double>(jobz, N, N, *h_A, *eigenvals,
        lapackf77_dsyevd);
  }

template<class T>
  int
  carma_syevd_cpu(char jobz, magma_int_t N, T *h_A, T *eigenvals) {
    MAGMA_TRACE("carma_getri not implemented with this type! \n");
    return EXIT_FAILURE;
  }
template<>
  int
  carma_syevd_cpu<float>(char jobz, magma_int_t N, float *h_A,
      float *eigenvals) {
    return carma_syevd_cpu_gen<float>(jobz, N, N, h_A, eigenvals,
        lapackf77_ssyevd);
  }
template<>
  int
  carma_syevd_cpu<double>(char jobz, magma_int_t N, double *h_A,
      double *eigenvals) {
    DEBUG_TRACE("ICI\n");
    return carma_syevd_cpu_gen<double>(jobz, N, N, h_A, eigenvals,
        lapackf77_dsyevd);
  }

template<class T>
  int
  carma_axpy_cpu(long N, T alpha, T *h_X, long incX, T *h_Y, long incY) {
    MAGMA_TRACE("carma_axpy_cpu not implemented with this type! \n");
    return EXIT_FAILURE;
  }
template<>
  int
  carma_axpy_cpu<float>(long N, float alpha, float *h_X, long incX, float *h_Y,
      long incY) {
    blasf77_saxpy(&N, &alpha, h_X, &incX, h_Y, &incY);
    return EXIT_SUCCESS;
  }
template<>
  int
  carma_axpy_cpu<double>(long N, double alpha, double *h_X, long incX,
      double *h_Y, long incY) {
    blasf77_daxpy(&N, &alpha, h_X, &incX, h_Y, &incY);
    return EXIT_SUCCESS;
  }

template<class T>
  int
  carma_gemm_cpu(char transa, char transb, long m, long n, long k, T alpha,
      T *A, long lda, T *B, long ldb, T beta, T *C, long ldc) {
    MAGMA_TRACE("carma_gemm_cpu not implemented with this type! \n");
    return EXIT_FAILURE;
  }
template<>
  int
  carma_gemm_cpu<float>(char transa, char transb, long m, long n, long k,
      float alpha, float *A, long lda, float *B, long ldb, float beta, float *C,
      long ldc) {
    blasf77_sgemm(&transa, &transb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta,
        C, &ldc);
    return EXIT_SUCCESS;
  }
template<>
  int
  carma_gemm_cpu<double>(char transa, char transb, long m, long n, long k,
      double alpha, double *A, long lda, double *B, long ldb, double beta,
      double *C, long ldc) {
    blasf77_dgemm(&transa, &transb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta,
        C, &ldc);
    return EXIT_SUCCESS;
  }

#else
#warning "MAGMA will not be used"
template<class T> int carma_svd(carma_obj<T> *imat, carma_obj<T> *eigenvals, carma_obj<T> *mod2act, carma_obj<T> *mes2mod)
{
  MAGMA_TRACE("!!!!!! MAGMA not compiled !!!!!!\n");
  return EXIT_FAILURE;
}
template int carma_svd<float>(carma_obj<float> *imat, carma_obj<float> *eigenvals, carma_obj<float> *mod2act, carma_obj<float> *mes2mod);
template int carma_svd<double>(carma_obj<double> *imat, carma_obj<double> *eigenvals, carma_obj<double> *mod2act, carma_obj<double> *mes2mod);

template<class T> int carma_syevd(char jobz, carma_obj<T> *mat, carma_host_obj<T> *eigenvals, carma_obj<T> *U)
{
  MAGMA_TRACE("!!!!!! MAGMA not compiled !!!!!!\n");
  return EXIT_FAILURE;
}
template int carma_syevd<float>(char jobz, carma_obj<float> *mat, carma_host_obj<float> *eigenvals, carma_obj<float> *U);
template int carma_syevd<double>(char jobz, carma_obj<double> *mat, carma_host_obj<double> *eigenvals, carma_obj<double> *U);

template<class T> int carma_syevd(char jobz, carma_obj<T> *mat,  carma_host_obj<T> *eigenvals)
{
  MAGMA_TRACE("!!!!!! MAGMA not compiled !!!!!!\n");
  return EXIT_FAILURE;
}
template int carma_syevd<float>(char jobz, carma_obj<float> *mat, carma_host_obj<float> *eigenvals);
template int carma_syevd<double>(char jobz, carma_obj<double> *mat, carma_host_obj<double> *eigenvals);

template<class T> int carma_syevd_m(long ngpu, char jobz, long N, T *mat, T *eigenvals)
{
  MAGMA_TRACE("!!!!!! MAGMA not compiled !!!!!!\n");
  return EXIT_FAILURE;
}
template int carma_syevd_m<float>(long ngpu, char jobz, long N, float *mat, float *eigenvals);
template int carma_syevd_m<double>(long ngpu, char jobz, long N, double *mat, double *eigenvals);

template<class T> int carma_getri(carma_obj<T> *d_iA)
{
  MAGMA_TRACE("!!!!!! MAGMA not compiled !!!!!!\n");
  return EXIT_FAILURE;
}
template int carma_getri<float>(carma_obj<float> *d_iA);
template int carma_getri<double>(carma_obj<double> *d_iA);

template<class T> int carma_potri(carma_obj<T> *d_iA)
{
  MAGMA_TRACE("!!!!!! MAGMA not compiled !!!!!!\n");
  return EXIT_FAILURE;
}
template int carma_potri<float>(carma_obj<float> *d_iA);
template int carma_potri<double>(carma_obj<double> *d_iA);

template<class T> int carma_potri(long num_gpus, T *h_A, carma_obj<T> *d_iA)
{
  MAGMA_TRACE("!!!!!! MAGMA not compiled !!!!!!\n");
  return EXIT_FAILURE;
}
template int carma_potri<float>(long num_gpus, float *h_A, carma_obj<float> *d_iA);
template int carma_potri<double>(long num_gpus, double *h_A, carma_obj<double> *d_iA);

template<class T_data> int carma_svd(carma_host_obj<T_data> *imat, carma_host_obj<T_data> *eigenvals, carma_host_obj<T_data> *mod2act, carma_host_obj<T_data> *mes2mod)
{
  MAGMA_TRACE("!!!!!! MAGMA not compiled !!!!!!\n");
  return EXIT_FAILURE;
}
template int carma_svd<float>(carma_host_obj<float> *imat, carma_host_obj<float> *eigenvals, carma_host_obj<float> *mod2act, carma_host_obj<float> *mes2mod);
template int carma_svd<double>(carma_host_obj<double> *imat, carma_host_obj<double> *eigenvals, carma_host_obj<double> *mod2act, carma_host_obj<double> *mes2mod);

template<class T> int carma_potri(carma_host_obj<T> *h_A) {
  MAGMA_TRACE("!!!!!! MAGMA not compiled !!!!!!\n");
  return EXIT_FAILURE;
}
template int carma_potri<float>(carma_host_obj<float> *h_A);
template int carma_potri<double>(carma_host_obj<double> *h_A);

template<class T> int carma_syevd(char jobz, carma_host_obj<T> *h_A, carma_host_obj<T> *eigenvals) {
  MAGMA_TRACE("!!!!!! MAGMA not compiled !!!!!!\n");
  return EXIT_FAILURE;
}
template int carma_syevd<float>(char jobz, carma_host_obj<float> *h_A, carma_host_obj<float> *eigenvals);
template int carma_syevd<double>(char jobz, carma_host_obj<double> *h_A, carma_host_obj<double> *eigenvals);

#endif
