#include <carma_obj.h>
#include <carma_host_obj.h>
#include "magma.h"
#include "magma_lapack.h"

#ifndef max
#define max(a,b)  (((a)<(b))?(b):(a))
#endif
#ifndef min
#define min(a,b)  (((a)>(b))?(b):(a))
#endif

/** These templates are used to select the proper Iamax executable from T_data*/
template<class T> int carma_gesvd(int m, int n, T *mat, T *eigenvals, T *U, T*VT, T*h_work, int lwork);
/**< Generic template for Iamax executable selection */
template<> int carma_gesvd<float>(int m, int n, float *mat, float *eigenvals, float *U, float *VT, float *h_work, int lwork) {
  magma_int_t info;
  magma_sgesvd('A', 'A', m, n, mat, m, eigenvals, U, m, VT, n, h_work, lwork, &info);
  //lapackf77_sgesvd("A", "A", &m, &n, mat, &m, eigenvals, U, &m, VT, &n, h_work, &lwork, &info);

  return info;
}
template<> int carma_gesvd<double>(int m, int n, double *mat, double *eigenvals, double *U, double *VT, double *h_work, int lwork) {
  magma_int_t info;
  magma_dgesvd('A', 'A', m, n, mat, m, eigenvals, U, m, VT, n, h_work, lwork, &info);
  return info;
}

template<class T> int carma_syevd(carma_obj<T> *mat, T *eigenvals, carma_obj<T> *U,
                                  magma_int_t (*ptr_syevd_gpu)(char, char, magma_int_t, T*, magma_int_t, T*, T*, magma_int_t, T*, magma_int_t, magma_int_t*, magma_int_t, magma_int_t*)) {
  magma_init();
  magma_int_t *iwork;
  magma_int_t N = mat->getDims(1);
  magma_int_t info, lwork, liwork, lda = N, aux_iwork[1];
  T aux_work[1];
  T *h_R, *h_work;

  //fprintf(stderr, "%s@%d : N=%d\n", __FILE__, __LINE__, N);
  magma_int_t ldda = N; //((N + 31) / 32) * 32;

  //fprintf(stderr, "%s@%d : ldda=%d\n", __FILE__, __LINE__, ldda);
  /* Query for workspace sizes */
  ptr_syevd_gpu( 'V', 'U',
          N, NULL, ldda, NULL,
          NULL, lda,
          aux_work,  -1,
          aux_iwork, -1,
          &info );
  if (info != 0)
    printf("%s@%d magma returned error %d: %s.\n", __FILE__, __LINE__, (int) info, magma_strerror(info));

  //fprintf(stderr, "%s@%d : I'm here\n", __FILE__, __LINE__);
  lwork  = (magma_int_t)aux_work[0];
  liwork = aux_iwork[0];
  posix_memalign( (void**) &iwork, 32, (liwork)*sizeof(int) );

  cudaMallocHost( (void**) &h_R, (N*lda)*sizeof(T) );
  cudaMallocHost( (void**) &h_work, (lwork)*sizeof(T) );

  //fprintf(stderr, "%s@%d : I'm here\n", __FILE__, __LINE__);
  /* Initialize the matrix */
  U->copy(mat, 1, 1);

  //fprintf(stderr, "%s@%d : I'm here\n", __FILE__, __LINE__);
  ptr_syevd_gpu('V', 'U',
      N, U->getData(), ldda, eigenvals,
      h_R, lda,
      h_work, lwork,
      iwork, liwork,
      &info);

  if (info != 0)
    printf("%s@%d magma returned error %d: %s.\n", __FILE__, __LINE__, (int) info, magma_strerror(info));
  //fprintf(stderr, "%s@%d : I'm here\n", __FILE__, __LINE__);

  cudaFreeHost(h_R);
  cudaFreeHost(h_work);
  free(iwork);
  magma_finalize();

  return EXIT_SUCCESS;
}

template<class T> int carma_syevd(carma_obj<T> *mat, T *eigenvals, carma_obj<T> *U);
template<> int carma_syevd<float>(caObjS *mat, float *eigenvals, caObjS *U){
  return carma_syevd<float>(mat, eigenvals, U, magma_ssyevd_gpu);
}
template<> int carma_syevd<double>(caObjD *mat, double *eigenvals, caObjD *U){
  return carma_syevd<double>(mat, eigenvals, U, magma_dsyevd_gpu);
}

template<class T, int method> int carma_syevd(carma_obj<T> *mat, T *eigenvals, carma_obj<T> *U);
template<> int carma_syevd<float,1>(caObjS *mat, float *eigenvals, caObjS *U){
  return carma_syevd<float>(mat, eigenvals, U, magma_ssyevd_gpu);
}
template<> int carma_syevd<double,1>(caObjD *mat, double *eigenvals, caObjD *U){
  return carma_syevd<double>(mat, eigenvals, U, magma_dsyevd_gpu);
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

template<class T> int carma_svd(carma_obj<T> *imat, carma_obj<T> *eigenvals, carma_obj<T> *mod2act, carma_obj<T> *mes2mod){
  //TODO: carma_svd
  cerr << "carma_svd not implemented ! \n";
 return 0;
}
template int carma_svd<float>(caObjS *imat, caObjS *eigenvals, caObjS *mod2act, caObjS *mes2mod);
template int carma_svd<double>(caObjD *imat, caObjD *eigenvals, caObjD *mod2act, caObjD *mes2mod);

template<class T> int carma_svd(carma_host_obj<T> *imat, carma_host_obj<T> *eigenvals, carma_host_obj<T> *mod2act, carma_host_obj<T> *mes2mod) {
  long int *dims = imat->getDims();
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

  carma_host_obj<T> *tmp = new carma_host_obj<T>(imat, MA_PAGELOCK);
  carma_gesvd(m, n, tmp->getData(), eigenvals->getData(), mes2mod->getData(), mod2act->getData(), h_work->getData(), lwork);

  delete h_work;
  delete tmp;

  return EXIT_SUCCESS;
}

template int carma_svd<float>(carma_host_obj<float> *imat, carma_host_obj<float> *eigenvals, carma_host_obj<float> *mod2act, carma_host_obj<float> *mes2mod);
template int carma_svd<double>(carma_host_obj<double> *imat, carma_host_obj<double> *eigenvals, carma_host_obj<double> *mod2act, carma_host_obj<double> *mes2mod);

