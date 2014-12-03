#ifndef _SUTRA_AO_UTILS_H_
#define _SUTRA_AO_UTILS_H_
#include <carma.h>
#include <carma_obj.h>

int
cfillrealp(cuFloatComplex *d_odata, float *d_idata, int N, int device);
int
cgetrealp(float *d_odata, cuFloatComplex *d_idata, int N, int device);
int
abs2(float *d_odata, cuFloatComplex *d_idata, int N, int device);
int
abs2c(cuFloatComplex *d_odata, cuFloatComplex *d_idata, int N, int device);
int
subap_norm(float *d_odata, float *d_idata, float *fact, float *norm,
    float nphot, int n, int N, int device);
int
fillindx(float *d_odata, float *d_idata, int *indx, float alpha, float beta,
    int N, int device);
int
fillindx(float *d_odata, float *d_idata, int *indx, float alpha, int N,
    int device);
int
fillindx(float *d_odata, float *d_idata, int *indx, int N, int device);
int
fillarr2d(float *d_odata, float *d_idata, int x0, int Ncol, int NC, int N,
    int device);
int
getarr2d(float *d_odata, float *d_idata, int x0, int Ncol, int NC, int N,
    int device);
template<class T>
int
addai(T *d_odata, T *i_data, int i, int sgn, int N, int device);
int
subap_norm_async(float *d_odata, float *d_idata, float *fact, float *norm,
    float nphot, int n, int N, carma_streams *streams, int device);
// templates
template<class T>
int
roll(T *idata, int N, int M, int nim);
template<class T>
int
roll(T *idata, int N, int M);
template<class T>
int
sutra_invgene(carma_obj<T> *imat, carma_obj<T> *cmat, carma_obj<T> *eigenvals,
    carma_obj<T> *mod2act, carma_obj<T> *mes2mod, int nfilt);
template<class T>
int
remove_avg(T *data, int N, int device);

#endif // _SUTRA_AO_UTILS_H_
