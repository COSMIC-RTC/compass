#include <yoga_obj.h>
#include <yoga_host_obj.h>
#include "magma.h"
#include "magma_lapack.h"

#ifndef max
#define max(a,b)  (((a)<(b))?(b):(a))
#endif
#ifndef min
#define min(a,b)  (((a)>(b))?(b):(a))
#endif

/** These templates are used to select the proper Iamax executable from T_data*/
template<class T> int yoga_sgesvd(int m, int n, T *mat, T *eigenvals, T *U, T*VT, T*h_work, int lwork);
/**< Generic template for Iamax executable selection */
template<> int yoga_sgesvd<float>(int m, int n, float *mat, float *eigenvals, float *U, float *VT, float *h_work, int lwork)
{
  magma_int_t  info;
  magma_sgesvd('A', 'A', m, n, mat, m, eigenvals, U, m, VT, n, h_work, lwork, &info);
  //lapackf77_sgesvd("A", "A", &m, &n, mat, &m, eigenvals, U, &m, VT, &n, h_work, &lwork, &info);

  return info;
}
template<> int yoga_sgesvd<double>(int m, int n, double *mat, double *eigenvals, double *U, double *VT, double *h_work, int lwork)
{
  magma_int_t info;
  magma_dgesvd('A', 'A', m, n, mat, m, eigenvals, U, m, VT, n, h_work, lwork, &info);
  return info;
}

template <class T> int yoga_svd(yoga_obj<T> *imat, yoga_obj<T> *eigenvals, yoga_obj<T> *mod2act, yoga_obj<T> *mes2mod)
{

	  int m = imat->getDims(1);
	  int n = imat->getDims(2);

  int n2  = m * n;
  int min_mn = m < n ? m : n;
  int max_mn = m > n ? m : n;
  
  long *dims_data = new long[2]; dims_data[0] = 1;
  dims_data[1] = n * n;
  yoga_host_obj<T> *VT    = new yoga_host_obj<T>(dims_data, MA_PAGELOCK);//MALLOC
  dims_data[1] = m * m;
  yoga_host_obj<T> *U    = new yoga_host_obj<T>(dims_data, MA_PAGELOCK);//MALLOC
  dims_data[1] = min_mn;
  yoga_host_obj<T> *S1    = new yoga_host_obj<T>(dims_data, MA_PAGELOCK);//MALLOC
  
  int nb    = 128; //magma_get_sgesvd_nb(n);
  int lwork = max(5*min_mn, (3*min_mn + max_mn))*nb;

  dims_data[1] = n2;
  yoga_host_obj<T> *h_A = new yoga_host_obj<T>(dims_data, MA_PAGELOCK);//PAGELOCK
  h_A->cpy_obj(imat, cudaMemcpyDeviceToHost);
  //imat->device2host(h_A->h_data);

  dims_data[1] = lwork;
  yoga_host_obj<T> *h_work = new yoga_host_obj<T>(dims_data, MA_PAGELOCK);//PAGELOCK
  //cutilSafeCall(cudaMallocHost( (void**)&h_work, lwork *sizeof(T) ));
  //cutilSafeCall(cudaMallocHost( (void**)&h_A, n2 *sizeof(T) ));

  yoga_sgesvd(m, n, h_A->getData(), S1->getData(), U->getData(), VT->getData(), h_work->getData(), lwork);
  
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

template int yoga_svd<float>(yObjS *imat, yObjS *eigenvals, yObjS *mod2act, yObjS *mes2mod);
template int yoga_svd<double>(yObjD *imat, yObjD *eigenvals, yObjD *mod2act, yObjD *mes2mod);


template <class T> int yoga_svd(yoga_host_obj<T> *imat, yoga_host_obj<T> *eigenvals, yoga_host_obj<T> *mod2act, yoga_host_obj<T> *mes2mod)
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
  yoga_host_obj<T> *h_work = new yoga_host_obj<T>(dims_data, MA_PAGELOCK);//PAGELOCK
  
  yoga_host_obj<T> *tmp = new yoga_host_obj<T>(imat, MA_PAGELOCK);
  yoga_sgesvd(m, n, tmp->getData(), eigenvals->getData(), mes2mod->getData(), mod2act->getData(), h_work->getData(), lwork);
  
  delete h_work;
  delete tmp;

  return EXIT_SUCCESS;
}

template int yoga_svd<float>(yoga_host_obj<float> *imat, yoga_host_obj<float> *eigenvals, yoga_host_obj<float> *mod2act, yoga_host_obj<float> *mes2mod);
template int yoga_svd<double>(yoga_host_obj<double> *imat, yoga_host_obj<double> *eigenvals, yoga_host_obj<double> *mod2act, yoga_host_obj<double> *mes2mod);

