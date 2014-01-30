#include <carma_obj.h>
#include <carma_host_obj.h>

#ifdef DEBUG
#define MAGMA_TRACE(fmt, args...) fprintf(stderr, fmt, ## args)
#else
#define MAGMA_TRACE(fmt, args...) /* */
#endif

#ifdef USE_MAGMA
#include "magma.h"
#include "magma_lapack.h"

#ifndef max
#define max(a,b)  (((a)<(b))?(b):(a))
#endif
#ifndef min
#define min(a,b)  (((a)>(b))?(b):(a))
#endif

template<class T> int carma_syevd(char jobz, carma_obj<T> *mat, T *eigenvals,
		magma_int_t (*ptr_syevd_gpu)(char jobz, char uplo, magma_int_t n, T *da,
				magma_int_t ldda, T *w, T *wa, magma_int_t ldwa, T *work,
				magma_int_t lwork, magma_int_t *iwork, magma_int_t liwork,
				magma_int_t *info)) {
	magma_int_t *iwork;
	magma_int_t N = mat->getDims(1);
	magma_int_t info, lwork, liwork, lda = N, aux_iwork[1];
	T aux_work[1];
	T *h_R, *h_work;

	//fprintf(stderr, "%s@%d : N=%d\n", __FILE__, __LINE__, N);
	magma_int_t ldda = N; //((N + 31) / 32) * 32;

	//fprintf(stderr, "%s@%d : ldda=%d\n", __FILE__, __LINE__, ldda);
	/* Query for workspace sizes */
	ptr_syevd_gpu(jobz, 'L', N, NULL, ldda, NULL, NULL, lda, aux_work, -1,
			aux_iwork, -1, &info);
	if (info != 0)
		printf("%s@%d magma returned error %d: %s.\n", __FILE__, __LINE__,
				(int) info, magma_strerror(info));

	//fprintf(stderr, "%s@%d : I'm here\n", __FILE__, __LINE__);
	lwork = (magma_int_t) aux_work[0];
	liwork = aux_iwork[0];
	magma_malloc_cpu((void**) &iwork, (liwork) * sizeof(int));

	cudaMallocHost((void**) &h_R, (N * lda) * sizeof(T));
	cudaMallocHost((void**) &h_work, (lwork) * sizeof(T));

	//fprintf(stderr, "%s@%d : I'm here\n", __FILE__, __LINE__);
	/* Initialize the matrix */

	//fprintf(stderr, "%s@%d : I'm here\n", __FILE__, __LINE__);
	ptr_syevd_gpu(jobz, 'L', N, mat->getData(), ldda, eigenvals, h_R, lda,
			h_work, lwork, iwork, liwork, &info);

	if (info != 0)
		printf("%s@%d magma returned error %d: %s.\n", __FILE__, __LINE__,
				(int) info, magma_strerror(info));
	//fprintf(stderr, "%s@%d : I'm here\n", __FILE__, __LINE__);

	cudaFreeHost(h_R);
	cudaFreeHost(h_work);
	free(iwork);

	return EXIT_SUCCESS;
}

template<class T> int carma_syevd(char jobz, carma_obj<T> *mat, T *eigenvals) {
	cerr << "carma_syevd: this method not implemented !\n";
	return EXIT_FAILURE;
}
template<> int carma_syevd<float>(char jobz, caObjS *mat, float *eigenvals) {
	return carma_syevd<float>(jobz, mat, eigenvals, magma_ssyevd_gpu);
}
template<> int carma_syevd<double>(char jobz, caObjD *mat, double *eigenvals) {
	return carma_syevd<double>(jobz, mat, eigenvals, magma_dsyevd_gpu);
}

template<class T> int carma_syevd(char jobz, carma_obj<T> *mat, T *eigenvals,
		carma_obj<T> *U) {
	cerr << "carma_syevd: this method not implemented !\n";
	return EXIT_FAILURE;
}
template<> int carma_syevd<float>(char jobz, caObjS *mat, float *eigenvals, caObjS *U) {
	U->copy(mat, 1, 1);
	return carma_syevd<float>(jobz, U, eigenvals, magma_ssyevd_gpu);
}
template<> int carma_syevd<double>(char jobz, caObjD *mat, double *eigenvals, caObjD *U) {
	U->copy(mat, 1, 1);
	return carma_syevd<double>(jobz, U, eigenvals, magma_dsyevd_gpu);
}

template<class T, int method> int carma_syevd(char jobz, carma_obj<T> *mat, T *eigenvals,
		carma_obj<T> *U) {
	cerr << "carma_syevd: this method not implemented !\n";
	return EXIT_FAILURE;
}
template<> int carma_syevd<float, 1>(char jobz, caObjS *mat, float *eigenvals, caObjS *U) {
	return carma_syevd<float>(jobz, mat, eigenvals, U);
}
template<> int carma_syevd<double, 1>(char jobz, caObjD *mat, double *eigenvals,
		caObjD *U) {
	return carma_syevd<double>(jobz, mat, eigenvals, U);
}
/*
 template<> int carma_syevd<float,2>(caObjS *mat, float *eigenvals, caObjS *U){
 return carma_syevd<float>(mat, eigenvals, U, magma_ssyevd_gpu_magmablas);
 }
 template<> int carma_syevd<double,2>(caObjD *mat, double *eigenvals, caObjD *U){
 return carma_syevd<double>(mat, eigenvals, U, magma_dsyevd_gpu_magmablas);
 }
 template<> int carma_syevd<float,3>(caObjS *mat, float *eigenvals, caObjS *U){
 return carma_syevd<float>(mat, eigenvals, U, magma_ssyevd_gpu_kblas);
 }
 template<> int carma_syevd<double,3>(caObjD *mat, double *eigenvals, caObjD *U){
 return carma_syevd<double>(mat, eigenvals, U, magma_dsyevd_gpu_kblas);
 }
 */

template<class T> int carma_syevd_m(magma_int_t ngpu, char jobz, T *mat, T *eigenvals,
		magma_int_t N,
		magma_int_t (*ptr_syevd_m)(magma_int_t nrgpu, char jobz, char uplo,
				magma_int_t n, T *a, magma_int_t lda, T *w, T *work,
				magma_int_t lwork, magma_int_t *iwork, magma_int_t liwork,
				magma_int_t *info)) {
	magma_int_t *iwork;
	magma_int_t info, lwork, liwork, lda = N, aux_iwork[1];
	T aux_work[1];
	T *h_R, *h_work;

	//fprintf(stderr, "%s@%d : ldda=%d\n", __FILE__, __LINE__, ldda);
	/* Query for workspace sizes */
	ptr_syevd_m(ngpu, jobz,  'L', N, NULL, lda, NULL, aux_work, -1, aux_iwork, -1,
			&info);
	if (info != 0)
		printf("%s@%d magma returned error %d: %s.\n", __FILE__, __LINE__,
				(int) info, magma_strerror(info));

	//fprintf(stderr, "%s@%d : I'm here\n", __FILE__, __LINE__);
	lwork = (magma_int_t) aux_work[0];
	liwork = aux_iwork[0];

	magma_malloc_cpu((void**) &iwork, (liwork) * sizeof(int));
	cudaMallocHost((void**) &h_R, (N * lda) * sizeof(T));
	cudaMallocHost((void**) &h_work, (lwork) * sizeof(T));

	//fprintf(stderr, "%s@%d : I'm here\n", __FILE__, __LINE__);
	ptr_syevd_m(ngpu, jobz,  'L', N, mat, lda, eigenvals, h_work, lwork, iwork,
			liwork, &info);
	//magma_dsyevd_m(opts.ngpu, opts.jobz, opts.uplo, N, h_R, lda, w1, h_work, lwork, iwork, liwork, &info );

	if (info != 0)
		printf("%s@%d magma returned error %d: %s.\n", __FILE__, __LINE__,
				(int) info, magma_strerror(info));
	//fprintf(stderr, "%s@%d : I'm here\n", __FILE__, __LINE__);

	cudaFreeHost(h_R);
	cudaFreeHost(h_work);
	free(iwork);

	return EXIT_SUCCESS;
}

template<class T> int carma_syevd_m(long ngpu, char jobz, T *mat, T *eigenvals, long N) {
	cerr << "carma_syevd_m not implemented this this type! \n";
	return EXIT_FAILURE;
}
template<> int carma_syevd_m<float>(long ngpu, char jobz, float *mat, float *eigenvals,
		long N) {
	return carma_syevd_m<float>(ngpu, jobz, mat, eigenvals, N, magma_ssyevd_m);
}
template<> int carma_syevd_m<double>(long ngpu, char jobz, double *mat, double *eigenvals,
		long N) {
	return carma_syevd_m<double>(ngpu, jobz, mat, eigenvals, N, magma_dsyevd_m);
}

template<class T> int carma_syevd_m(long ngpu, char jobz, T *mat, T *eigenvals, T *U,
		long N) {
	cerr << "carma_syevd_m not implemented this this type! \n";
	return EXIT_FAILURE;
}
template<> int carma_syevd_m<float>(long ngpu, char jobz, float *mat, float *eigenvals,
		float *U, long N) {
	/* Initialize the matrix */
	magma_int_t lda = N, n = N;
	lapackf77_slacpy(MagmaUpperLowerStr, &n, &n, mat, &lda, U, &lda);

	return carma_syevd_m<float>(ngpu, jobz, U, eigenvals, N, magma_ssyevd_m);
}
template<> int carma_syevd_m<double>(long ngpu, char jobz, double *mat, double *eigenvals,
		double *U, long N) {
	/* Initialize the matrix */
	magma_int_t lda = N, n = N;
	lapackf77_dlacpy(MagmaUpperLowerStr, &n, &n, mat, &lda, U, &lda);

	return carma_syevd_m<double>(ngpu, jobz, U, eigenvals, N, magma_dsyevd_m);
}

template<class T> int carma_svd(carma_obj<T> *mat, carma_obj<T> *eigenvals,
		carma_obj<T> *mod2act, carma_obj<T> *mes2mod) {
	//TODO: carma_svd
	cerr << "carma_svd not implemented this this type! \n";
	return EXIT_FAILURE;
}

template<> int carma_svd<float>(caObjS *imat, caObjS *eigenvals,
		caObjS *mod2act, caObjS *mes2mod) {
	cerr << "carma_svd not implemented on device object! \n";
	return EXIT_FAILURE;
	//return carma_gesvd<float>(mat, eigenvals, U, magma_sgesvd);
}
template<> int carma_svd<double>(caObjD *imat, caObjD *eigenvals,
		caObjD *mod2act, caObjD *mes2mod) {
	cerr << "carma_svd not implemented on device object! \n";
	return EXIT_FAILURE;
	//return carma_gesvd<double>(mat, eigenvals, U, magma_dgesvd);
}

template<class T> int carma_svd(carma_host_obj<T> *mat,
		carma_host_obj<T> *eigenvals, carma_host_obj<T> *mod2act,
		carma_host_obj<T> *mes2mod,
		magma_int_t (*ptr_gesvd)(char jobu, char jobvt, magma_int_t m,
				magma_int_t n, T *A, magma_int_t lda, T *s, T *U,
				magma_int_t ldu, T *VT, magma_int_t ldvt, T *work,
				magma_int_t lwork, magma_int_t *info)) {
	long int *dims = mat->getDims();
	int m = dims[1];
	int n = dims[2];
	int min_mn = m < n ? m : n;
	int max_mn = m > n ? m : n;

	long *dims_data = new long[2];
	dims_data[0] = 1;
	int nb = 128; //magma_get_sgesvd_nb(n);
	int lwork = max(5*min_mn, (3*min_mn + max_mn)) * nb;
	dims_data[1] = lwork;
	carma_host_obj<T> *h_work = new carma_host_obj<T>(dims_data, MA_PAGELOCK); //PAGELOCK

	carma_host_obj<T> *tmp = new carma_host_obj<T>(mat, MA_PAGELOCK);
	//carma_gesvd(m, n, tmp->getData(), eigenvals->getData(), mes2mod->getData(), mod2act->getData(), h_work->getData(), lwork);
	magma_int_t info;
	ptr_gesvd('A', 'A', m, n, tmp->getData(), m, eigenvals->getData(),
			mes2mod->getData(), m, mod2act->getData(), n, h_work->getData(),
			lwork, &info);

	delete h_work;
	delete tmp;

	return EXIT_SUCCESS;
}

template<class T> int carma_svd(carma_host_obj<T> *mat,
		carma_host_obj<T> *eigenvals, carma_host_obj<T> *mod2act,
		carma_host_obj<T> *mes2mod) {
	//TODO: carma_svd
	cerr << "carma_svd not implemented with this type! \n";
	return EXIT_FAILURE;
}
template<> int carma_svd<float>(carma_host_obj<float> *mat,
		carma_host_obj<float> *eigenvals, carma_host_obj<float> *mod2act,
		carma_host_obj<float> *mes2mod) {
	return carma_svd(mat, eigenvals, mod2act, mes2mod, magma_sgesvd);
}
template<> int carma_svd<double>(carma_host_obj<double> *mat,
		carma_host_obj<double> *eigenvals, carma_host_obj<double> *mod2act,
		carma_host_obj<double> *mes2mod) {
	return carma_svd(mat, eigenvals, mod2act, mes2mod, magma_dgesvd);
}

template<class T> int carma_potri(carma_obj<T> *d_iA,
		magma_int_t (*ptr_potrf)(char uplo, magma_int_t n, T *d_A,
				magma_int_t ldda, magma_int_t *info),
		magma_int_t (*ptr_potri)(char uplo, magma_int_t n, T *d_A,
				magma_int_t ldda, magma_int_t *info)) {
	long int *dims = d_iA->getDims();
	int m = dims[1];
	int n = dims[2];
	if (m != n) {
		cerr << "carma_potri : non symmetric matrix\n";
		return EXIT_FAILURE;
	}

	magma_int_t info;
	ptr_potrf('U', n, d_iA->getData(), n, &info);
	ptr_potri('U', n, d_iA->getData(), n, &info);
	//ptr_potrf('L', n, d_iA->getData(), n, &info);
	//ptr_potri('L', n, d_iA->getData(), n, &info);

	fill_sym_matrix('L', d_iA->getData(), n, n * m);
	return EXIT_SUCCESS;
}

template<class T> int carma_potri(carma_obj<T> *d_iA) {
	cerr << "carma_potri not implemented with this type! \n";
	return EXIT_FAILURE;
}
template<> int carma_potri<float>(carma_obj<float> *d_iA) {
	return carma_potri(d_iA, magma_spotrf_gpu, magma_spotri_gpu);
}
template<> int carma_potri<double>(carma_obj<double> *d_iA) {
	return carma_potri(d_iA, magma_dpotrf_gpu, magma_dpotri_gpu);
}

template<class T> int carma_potri(long num_gpus, T *h_A, carma_obj<T> *d_iA,
		magma_int_t (*ptr_potrf)(magma_int_t num_gpus, char uplo, magma_int_t n,
				T **d_A, magma_int_t ldda, magma_int_t *info),
		magma_int_t (*ptr_potri)(char uplo, magma_int_t n, T *d_A,
				magma_int_t ldda, magma_int_t *info),
		magma_int_t (*ptr_get_potrf_nb)(magma_int_t m),
		void (*ptr_setmatrix_1D_col_bcyclic)(magma_int_t m, magma_int_t n,
				const T *hA, magma_int_t lda, T *dA[], magma_int_t ldda,
				magma_int_t num_gpus, magma_int_t nb),
		void (*ptr_getmatrix_1D_col_bcyclic)(magma_int_t m, magma_int_t n,
				T *dA[], magma_int_t ldda, T *hA, magma_int_t lda,
				magma_int_t num_gpus, magma_int_t nb)) {
	long int *dims = d_iA->getDims();
	int m = dims[1];
	int n = dims[2];
	if (m != n) {
		cerr << "carma_potri : non symmetric matrix\n";
		return EXIT_FAILURE;
	}
	magma_int_t nb = ptr_get_potrf_nb(n);
	magma_int_t ldda = ((n + nb - 1) / nb) * nb;
	magma_int_t lda = n;

	T *d_lA[num_gpus];
	magma_int_t max_size = nb * (1 + n / (nb * num_gpus)) * nb
			* ((n + nb - 1) / nb);
	for (int dev = 0; dev < num_gpus; dev++) {
		magma_setdevice(dev);
		cudaMalloc(&d_lA[dev], max_size * sizeof(T));
	}
	ptr_setmatrix_1D_col_bcyclic(n, n, h_A, lda, d_lA, ldda, num_gpus, nb);

	magma_int_t info;
	ptr_potrf(num_gpus, 'U', n, d_lA, n, &info);
	ptr_getmatrix_1D_col_bcyclic(n, n, d_lA, ldda, h_A, lda, num_gpus, nb);

	d_iA->host2device(h_A);
	ptr_potri('U', n, d_iA->getData(), n, &info);

	fill_sym_matrix('L', d_iA->getData(), n, n * m);
	return EXIT_SUCCESS;
}

template<class T> int carma_potri(long num_gpus, T *h_A, carma_obj<T> *d_iA) {
	cerr << "carma_potri (multiGPU) not implemented with this type! \n";
	return EXIT_FAILURE;
}
/*
 template<> int carma_potri<float>(int num_gpus, float *h_A, carma_obj<float> *d_iA) {
 return carma_potri(num_gpus, h_A, d_iA, magma_spotrf_mgpu, magma_spotri_gpu, magma_get_spotrf_nb, magma_ssetmatrix_1D_col_bcyclic, magma_sgetmatrix_1D_col_bcyclic);
 }
 template<> int carma_potri<double>(int num_gpus, double *h_A, carma_obj<double> *d_iA) {
 return carma_potri(num_gpus, h_A, d_iA, magma_dpotrf_mgpu, magma_dpotri_gpu, magma_get_dpotrf_nb, magma_dsetmatrix_1D_col_bcyclic, magma_dgetmatrix_1D_col_bcyclic);
 }
 */
template<class T> int carma_getri(carma_obj<T> *d_iA,
		magma_int_t (*ptr_getrf)(magma_int_t m, magma_int_t n, T *dA,
				magma_int_t ldda, magma_int_t *ipiv, magma_int_t *info),
		magma_int_t (*ptr_getri)(magma_int_t n, T *dA, magma_int_t ldda,
				magma_int_t *ipiv, T *dwork, magma_int_t lwork,
				magma_int_t *info),
		magma_int_t (*ptr_get_getri_nb)(magma_int_t n)) {
	long int *dims = d_iA->getDims();
	int m = dims[1];
	int n = dims[2];

	if (m != n) {
		cerr << "carma_getri : non symmetric matrix\n";
		return EXIT_FAILURE;
	}
	magma_int_t *ipiv;
	magma_malloc_cpu((void**) &ipiv, n * sizeof(int));

	T *dwork;
	magma_int_t lda = n;
	magma_int_t ldda = ((n + 31) / 32) * 32;
	magma_int_t ldwork = n * ptr_get_getri_nb(n);
	if (lda * ldda != d_iA->getNbElem())
		cerr << "Please provide arrays of proper size\n";

	magma_int_t info;
	cudaMalloc((void**) &dwork, ldwork * sizeof(T));
	ptr_getrf(n, n, d_iA->getData(), ldda, ipiv, &info);
	ptr_getri(n, d_iA->getData(), ldda, ipiv, dwork, ldwork, &info);

	cudaFree(dwork);
	free(ipiv);

	return EXIT_SUCCESS;
}

template<class T> int carma_getri(carma_obj<T> *d_iA) {
	cerr << "carma_getri not implemented with this type! \n";
	return EXIT_FAILURE;
}
template<> int carma_getri<float>(carma_obj<float> *d_iA) {
	return carma_getri(d_iA, magma_sgetrf_gpu, magma_sgetri_gpu,
			magma_get_sgetri_nb);
}
template<> int carma_getri<double>(carma_obj<double> *d_iA) {
	return carma_getri(d_iA, magma_dgetrf_gpu, magma_dgetri_gpu,
			magma_get_dgetri_nb);
}

template<class T> int carma_getri(T* h_A, carma_obj<T> *d_iA) {
	cerr << "carma_potrfsi not implemented with this type! \n";
	return EXIT_FAILURE;
}
template<> int carma_getri<float>(float* h_A, carma_obj<float> *d_iA) {
	long int *dims = d_iA->getDims();
	int N = dims[2];
	magma_int_t lda = N;
	magma_int_t ldda = ((N + 31) / 32) * 32;
	magma_ssetmatrix(N, N, h_A, lda, d_iA->getData(), ldda);

	return carma_getri(d_iA, magma_sgetrf_gpu, magma_sgetri_gpu,
			magma_get_sgetri_nb);
}
template<> int carma_getri<double>(double* h_A, carma_obj<double> *d_iA) {
	long int *dims = d_iA->getDims();
	int N = dims[2];
	magma_int_t lda = N;
	magma_int_t ldda = ((N + 31) / 32) * 32;
	magma_dsetmatrix(N, N, h_A, lda, d_iA->getData(), ldda);

	return carma_getri(d_iA, magma_dgetrf_gpu, magma_dgetri_gpu,
			magma_get_dgetri_nb);
}

#else
#warning "MAGMA will not be used"
template<class T> int carma_svd(carma_obj<T> *imat, carma_obj<T> *eigenvals, carma_obj<T> *mod2act, carma_obj<T> *mes2mod) {
  MAGMA_TRACE("!!!!!! MAGMA not compiled !!!!!!\n");
  return EXIT_FAILURE;
}
template int carma_svd<float>(carma_obj<float> *imat, carma_obj<float> *eigenvals, carma_obj<float> *mod2act, carma_obj<float> *mes2mod);
template int carma_svd<double>(carma_obj<double> *imat, carma_obj<double> *eigenvals, carma_obj<double> *mod2act, carma_obj<double> *mes2mod);

template<class T> int carma_syevd(char jobz, carma_obj<T> *mat, T *eigenvals, carma_obj<T> *U){
  MAGMA_TRACE("!!!!!! MAGMA not compiled !!!!!!\n");
  return EXIT_FAILURE;
}
template int carma_syevd<float>(char jobz, carma_obj<float> *mat, float *eigenvals, carma_obj<float> *U);
template int carma_syevd<double>(char jobz, carma_obj<double> *mat, double *eigenvals, carma_obj<double> *U);

template<class T> int carma_syevd(char jobz, carma_obj<T> *mat, T *eigenvals){
  MAGMA_TRACE("!!!!!! MAGMA not compiled !!!!!!\n");
  return EXIT_FAILURE;
}
template int carma_syevd<float>(char jobz, carma_obj<float> *mat, float *eigenvals);
template int carma_syevd<double>(char jobz, carma_obj<double> *mat, double *eigenvals);

template<class T> int carma_syevd_m(int ngpu, char jobz, T *mat, T *eigenvals, T *U, int N){
  MAGMA_TRACE("!!!!!! MAGMA not compiled !!!!!!\n");
  return EXIT_FAILURE;
}
template int carma_syevd_m<float>(int ngpu, char jobz, float *mat, float *eigenvals, float *U, int N);
template int carma_syevd_m<double>(int ngpu, char jobz, double *mat, double *eigenvals, double *U, int N);

template<class T> int carma_getri(carma_obj<T> *d_iA) {
  MAGMA_TRACE("!!!!!! MAGMA not compiled !!!!!!\n");
  return EXIT_FAILURE;
}
template int carma_getri<float>(carma_obj<float> *d_iA);
template int carma_getri<double>(carma_obj<double> *d_iA);

template<class T> int carma_getri(T *h_A, carma_obj<T> *d_iA) {
  MAGMA_TRACE("!!!!!! MAGMA not compiled !!!!!!\n");
  return EXIT_FAILURE;
}
template int carma_getri<float>(float *h_A, carma_obj<float> *d_iA);
template int carma_getri<double>(double *h_A, carma_obj<double> *d_iA);

template<class T> int carma_potri(carma_obj<T> *d_iA) {
  MAGMA_TRACE("!!!!!! MAGMA not compiled !!!!!!\n");
  return EXIT_FAILURE;
}
template int carma_potri<float>(carma_obj<float> *d_iA);
template int carma_potri<double>(carma_obj<double> *d_iA);

template<class T> int carma_potri(long num_gpus, T *h_A, carma_obj<T> *d_iA){
  MAGMA_TRACE("!!!!!! MAGMA not compiled !!!!!!\n");
  return EXIT_FAILURE;
}
template int carma_potri<float>(long num_gpus, float *h_A, carma_obj<float> *d_iA);
template int carma_potri<double>(long num_gpus, double *h_A, carma_obj<double> *d_iA);

template<class T_data> int carma_svd(carma_host_obj<T_data> *imat, carma_host_obj<T_data> *eigenvals, carma_host_obj<T_data> *mod2act, carma_host_obj<T_data> *mes2mod) {
  MAGMA_TRACE("!!!!!! MAGMA not compiled !!!!!!\n");
  return EXIT_FAILURE;
}
template int carma_svd<float>(carma_host_obj<float> *imat, carma_host_obj<float> *eigenvals, carma_host_obj<float> *mod2act, carma_host_obj<float> *mes2mod);
template int carma_svd<double>(carma_host_obj<double> *imat, carma_host_obj<double> *eigenvals, carma_host_obj<double> *mod2act, carma_host_obj<double> *mes2mod);
#endif
