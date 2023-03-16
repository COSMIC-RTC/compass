// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_groot.h
//! \ingroup   libsutra
//! \class     SutraGroot
//! \brief     this class provides the groot features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.0
//! \date      2022/01/24

#ifndef _SUTRA_GROOT_H_
#define _SUTRA_GROOT_H_

#include <carma.h>
#include <carma_cublas.h>
#include <carma_host_obj.h>
#include <carma_obj.h>
#include <carma_sparse_obj.h>

class SutraGroot {
 public:
  CarmaContext *current_context;
  int device;

  int nactus;     // number of actuators
  int nssp;       // number of valid ssp
  int nlayers;    // number of atmos layers
  int npts;       // number of samples for aliasig computation
  float gsangle;  // Guide star angle [rad]
  float fc;       // DM cut-off frequency [m]
  float scale;
  float d;

  CarmaHostObj<float> *h_vdt;      // v*dt/g [m/s]
  CarmaHostObj<float> *h_Htheta;   // H*theta (theta = GS radial distance) [m]
  CarmaHostObj<float> *h_L0;       // L0 [m]
  CarmaHostObj<float> *h_winddir;  // wind directions [rad]
  CarmaHostObj<float> *h_scale;    // r0**(-5/3) * frac * (lambda/2pi)**2
  CarmaHostObj<float> *h_weights;  // Simpson method weights for integration
  // Covariance matrix estimation
  CarmaObj<float> *d_Cerr;  // Residual error covariance matrix on DM actuators
  CarmaObj<float>
      *d_TT;  // TT component of the residual error covariance matrix
  CarmaObj<float> *d_CaXX;  // XX aliasing covariance matrix
  CarmaObj<float> *d_CaYY;  // YY aliasing covariance matrix

  CarmaObj<float>
      *d_TTPfilter;  // Tip-tilt and piston filter matrix (= Btt.dot(P))
  CarmaObj<float> *d_pzt2tt;  // pzt to TT matrix
  CarmaObj<float> *d_Nact;    // Coupling matrix
  CarmaObj<float> *d_xpos;    // X-positions of DM actuators or ssp [m]
  CarmaObj<float> *d_ypos;    // Y-positions of DM actuators or ssp [m]

  // Dphi lowpass
  CarmaObj<float> *d_tab_int_x;  // Tabulated integral
  CarmaObj<float> *d_tab_int_y;

 public:
  SutraGroot(CarmaContext *context, int device, int nactus, int nlayers,
              float gsangle, float *vdt, float *Htheta, float *L0,
              float *winddir, float *scale, float *pzt2tt, float *TTPfilter,
              float *Nact, float *xpos, float *ypos, float fc);

  SutraGroot(CarmaContext *context, int device, int nssp, float *weights,
              float scale, float *xpos, float *ypos, float fc, float d,
              int npts);
  ~SutraGroot();

  cublasHandle_t cublas_handle() { return current_context->get_cublas_handle(); }
  void init_common(CarmaContext *context, int device, float *xpos, float *ypos,
                   int N, float fc);
  int compute_Cerr();
  int compute_Calias();
};
template <class T_data>
int compute_Cerr_layer(T_data *Cerr, int N, T_data *tab_int_x,
                       T_data *tab_int_y, T_data *xpos, T_data *ypos,
                       T_data vdt, T_data Htheta, T_data L0, T_data fc,
                       T_data winddir, T_data gsangle, T_data scale, int Ntab,
                       CarmaDevice *device);
template <class T_data>
int compute_Ca(T_data *CaXX, T_data *CaYY, int nssp, T_data *tab_int_x,
               T_data *tab_int_y, T_data *xpos, T_data *ypos, T_data offset,
               T_data d, T_data fc, T_data scale, T_data weight, int Ntab,
               CarmaDevice *device);
template <class T_data>
int tab_u831J0(T_data *tab_int_x, T_data *tab_int_y, int N,
               CarmaDevice *device);
template <class T_data>
void cumsum(T_data *odata, T_data *idata, int N);
template <class T_data>
int add_transpose(T_data *Cerr, int N, CarmaDevice *device);

#endif
