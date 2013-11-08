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
template<class T> int carma_sgesvd(int m, int n, T *mat, T *eigenvals, T *U, T*VT, T*h_work, int lwork);
/**< Generic template for Iamax executable selection */
template<> int carma_sgesvd<float>(int m, int n, float *mat, float *eigenvals, float *U, float *VT, float *h_work, int lwork)
{
  magma_int_t  info;
  magma_sgesvd('A', 'A', m, n, mat, m, eigenvals, U, m, VT, n, h_work, lwork, &info);
  //lapackf77_sgesvd("A", "A", &m, &n, mat, &m, eigenvals, U, &m, VT, &n, h_work, &lwork, &info);

  return info;
}
template<> int carma_sgesvd<double>(int m, int n, double *mat, double *eigenvals, double *U, double *VT, double *h_work, int lwork)
{
  magma_int_t info;
  magma_dgesvd('A', 'A', m, n, mat, m, eigenvals, U, m, VT, n, h_work, lwork, &info);
  return info;
}

template <class T> int carma_svd(carma_obj<T> *imat, carma_obj<T> *eigenvals, carma_obj<T> *mod2act, carma_obj<T> *mes2mod)
{

	  int m = imat->getDims(1);
	  int n = imat->getDims(2);

  int n2  = m * n;
  int min_mn = m < n ? m : n;
  int max_mn = m > n ? m : n;
  
  long *dims_data = new long[2]; dims_data[0] = 1;
  dims_data[1] = n * n;
  carma_host_obj<T> *VT    = new carma_host_obj<T>(dims_data, MA_PAGELOCK);//MALLOC
  dims_data[1] = m * m;
  carma_host_obj<T> *U    = new carma_host_obj<T>(dims_data, MA_PAGELOCK);//MALLOC
  dims_data[1] = min_mn;
  carma_host_obj<T> *S1    = new carma_host_obj<T>(dims_data, MA_PAGELOCK);//MALLOC
  
  int nb    = 128; //magma_get_sgesvd_nb(n);
  int lwork = max(5*min_mn, (3*min_mn + max_mn))*nb;

  dims_data[1] = n2;
  carma_host_obj<T> *h_A = new carma_host_obj<T>(dims_data, MA_PAGELOCK);//PAGELOCK
  h_A->cpy_obj(imat, cudaMemcpyDeviceToHost);
  //imat->device2host(h_A->h_data);

  dims_data[1] = lwork;
  carma_host_obj<T> *h_work = new carma_host_obj<T>(dims_data, MA_PAGELOCK);//PAGELOCK
  //cutilSafeCall(cudaMallocHost( (void**)&h_work, lwork *sizeof(T) ));
  //cutilSafeCall(cudaMallocHost( (void**)&h_A, n2 *sizeof(T) ));

  carma_sgesvd(m, n, h_A->getData(), S1->getData(), U->getData(), VT->getData(), h_work->getData(), lwork);
  
  VT->cpy_obj(mod2act, cudaMemcpyHostToDevice);
  U->cpy_obj(mes2mod, cudaMemcpyHostToDevice);
  S1->cpy_obj(eigenvals, cudaMemcpyHostToDevice);
//  mod2act->host2device(VT->h_data);
//  mes2mod->host2device(U->h_data);
//  eigenvals->host2device(S1->h_data);
  
  //free(h_A);
  delete U;
  delete VT;
  delete S1;
  delete h_work;
  delete h_A;

  return EXIT_SUCCESS;
}

template int carma_svd<float>(caObjS *imat, caObjS *eigenvals, caObjS *mod2act, caObjS *mes2mod);
template int carma_svd<double>(caObjD *imat, caObjD *eigenvals, caObjD *mod2act, caObjD *mes2mod);


template <class T> int carma_svd(carma_host_obj<T> *imat, carma_host_obj<T> *eigenvals, carma_host_obj<T> *mod2act, carma_host_obj<T> *mes2mod)
{
	long int *dims = imat->getDims();
  int m = dims[1];
  int n = dims[2];
  int min_mn = m < n ? m : n;
  int max_mn = m > n ? m : n;
  
  long *dims_data = new long[2]; dims_data[0] = 1;
  int nb    = 128; //magma_get_sgesvd_nb(n);
  int lwork = max(5*min_mn, (3*min_mn + max_mn))*nb;
  dims_data[1] = lwork;
  carma_host_obj<T> *h_work = new carma_host_obj<T>(dims_data, MA_PAGELOCK);//PAGELOCK
  
  carma_host_obj<T> *tmp = new carma_host_obj<T>(imat, MA_PAGELOCK);
  carma_sgesvd(m, n, tmp->getData(), eigenvals->getData(), mes2mod->getData(), mod2act->getData(), h_work->getData(), lwork);
  
  delete h_work;
  delete tmp;

  return EXIT_SUCCESS;
}

template int carma_svd<float>(carma_host_obj<float> *imat, carma_host_obj<float> *eigenvals, carma_host_obj<float> *mod2act, carma_host_obj<float> *mes2mod);
template int carma_svd<double>(carma_host_obj<double> *imat, carma_host_obj<double> *eigenvals, carma_host_obj<double> *mod2act, carma_host_obj<double> *mes2mod);

