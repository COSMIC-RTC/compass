/**
 * \file sutra_utils.h
 *
 * \ingroup libsutra
 *
 * \brief this class provides utilities to COMPASS
 *
 * \authors Damien Gratadour & Arnaud Sevin & Florian Ferreira
 *
 * \version 1.0
 *
 * \date 2011/01/28
 *
 */
#ifndef _SUTRA_UTILS_H_
#define _SUTRA_UTILS_H_
#include <carma.h>
#include <carma_obj.h>
#include <carma_sparse_obj.h>

int compute_nmaxhr(long nvalid);
int cfillrealp(cuFloatComplex *d_odata, float *d_idata, int N,
               carma_device *device);
int cgetrealp(float *d_odata, cuFloatComplex *d_idata, int N,
              carma_device *device);
int abs2(float *d_odata, cuFloatComplex *d_idata, int N, carma_device *device);
int abs2(float *d_odata, cuFloatComplex *d_idata, int N, float fact,
         carma_device *device);
int abs2c(cuFloatComplex *d_odata, cuFloatComplex *d_idata, int N,
          carma_device *device);
int convolve(cuFloatComplex *d_odata, cuFloatComplex *d_idata, int N,
             carma_device *device);
int convolve_modulate(cuFloatComplex *d_odata, cuFloatComplex *d_idata, int mod,
                      int N, carma_device *device);
int subap_norm(float *d_odata, float *d_idata, float *fact, float *norm,
               float nphot, int n, int N, carma_device *device);
int fillindx(float *d_odata, float *d_idata, int *indx, float alpha, float beta,
             int N, carma_device *device);
int fillindx(float *d_odata, float *d_idata, int *indx, float alpha, int N,
             carma_device *device);
int fillindx(float *d_odata, float *d_idata, int *indx, int N,
             carma_device *device);
int fillarr2d(float *d_odata, float *d_idata, int x0, int Ncol, int NC, int N,
              carma_device *device);
int fillarr2d(float *d_odata, float *d_idata, int x0, int Ncol, int NC, int N,
              int dir, carma_device *device);
int getarr2d(float *d_odata, float *d_idata, int x0, int Ncol, int NC, int N,
             carma_device *device);
template <class T>
int addai(T *d_odata, T *i_data, int i, int sgn, int N, carma_device *device);
int subap_norm_async(float *d_odata, float *d_idata, float *fact, float *norm,
                     float nphot, int n, int N, carma_streams *streams,
                     carma_device *device);
// templates
template <class T>
int roll(T *idata, int N, int M, int nim, carma_device *device);
template <class T>
int roll(T *idata, int N, int M, carma_device *device);
template <class T>
int roll_mult(T *odata, T *idata, int N, int M, T alpha, carma_device *device);
template <class T>
int sutra_invgene(carma_obj<T> *imat, carma_obj<T> *cmat,
                  carma_obj<T> *eigenvals, carma_obj<T> *mod2act,
                  carma_obj<T> *mes2mod, int nfilt);
template <class T>
int remove_avg(T *data, int N, carma_device *device);
template <class T>
int mult_vect(T *d_data, T *scale, int N, carma_device *device);
template <class T>
int mult_vect(T *d_data, T *scale, T gain, int N, carma_device *device);
template <class T>
int mult_vect(T *d_data, T gain, int N, carma_device *device);

int mult_int(float *o_data, float *i_data, float *scale, float gain, int N,
             carma_device *device);
int mult_int(float *o_data, float *i_data, float *scale, float gain, int N,
             carma_device *device, carma_streams *streams);
int mult_int(float *o_data, float *i_data, float gain, int N,
             carma_device *device);
int add_md(float *o_matrix, float *i_matrix, float *i_vector, int N,
           carma_device *device);

#endif  // _SUTRA_UTILS_H_
