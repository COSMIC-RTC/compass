// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_groot.hpp
//! \ingroup   libsutra
//! \class     SutraGroot
//! \brief     this class provides the groot features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#ifndef _SUTRA_GROOT_H_
#define _SUTRA_GROOT_H_

#include <carma.hpp>
#include <carma_cublas.hpp>
#include <carma_host_obj.hpp>
#include <carma_obj.hpp>
#include <carma_sparse_obj.hpp>

class SutraGroot {
 public:
  CarmaContext *current_context;
  int32_t device;

  int32_t nactus;     // number of actuators
  int32_t nssp;       // number of valid ssp
  int32_t nlayers;    // number of atmos layers
  int32_t npts;       // number of samples for aliasig computation
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
  SutraGroot(CarmaContext *context, int32_t device, int32_t nactus, int32_t nlayers,
              float gsangle, float *vdt, float *Htheta, float *L0,
              float *winddir, float *scale, float *pzt2tt, float *TTPfilter,
              float *Nact, float *xpos, float *ypos, float fc);

  SutraGroot(CarmaContext *context, int32_t device, int32_t nssp, float *weights,
              float scale, float *xpos, float *ypos, float fc, float d,
              int32_t npts);
  ~SutraGroot();

  cublasHandle_t cublas_handle() { return current_context->get_cublas_handle(); }
  void init_common(CarmaContext *context, int32_t device, float *xpos, float *ypos,
                   int32_t N, float fc);
  int32_t compute_Cerr();
  int32_t compute_Calias();
};
template <class T_data>
int32_t compute_Cerr_layer(T_data *Cerr, int32_t N, T_data *tab_int_x,
                       T_data *tab_int_y, T_data *xpos, T_data *ypos,
                       T_data vdt, T_data Htheta, T_data L0, T_data fc,
                       T_data winddir, T_data gsangle, T_data scale, int32_t Ntab,
                       CarmaDevice *device);
template <class T_data>
int32_t compute_Ca(T_data *CaXX, T_data *CaYY, int32_t nssp, T_data *tab_int_x,
               T_data *tab_int_y, T_data *xpos, T_data *ypos, T_data offset,
               T_data d, T_data fc, T_data scale, T_data weight, int32_t Ntab,
               CarmaDevice *device);
template <class T_data>
int32_t tab_u831J0(T_data *tab_int_x, T_data *tab_int_y, int32_t N,
               CarmaDevice *device);
template <class T_data>
void cumsum(T_data *odata, T_data *idata, int32_t N);
template <class T_data>
int32_t add_transpose(T_data *Cerr, int32_t N, CarmaDevice *device);

#endif
