// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2019 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.
//  Distributed under GNU - LGPL
//
//  COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
//  General Public License as published by the Free Software Foundation, either version 3 of the License,
//  or any later version.
//
//  COMPASS: End-to-end AO simulation tool using GPU acceleration
//  The COMPASS platform was designed to meet the need of high-performance for the simulation of AO systems.
//
//  The final product includes a software package for simulating all the critical subcomponents of AO,
//  particularly in the context of the ELT and a real-time core based on several control approaches,
//  with performances consistent with its integration into an instrument. Taking advantage of the specific
//  hardware architecture of the GPU, the COMPASS tool allows to achieve adequate execution speeds to
//  conduct large simulation campaigns called to the ELT.
//
//  The COMPASS platform can be used to carry a wide variety of simulations to both testspecific components
//  of AO of the E-ELT (such as wavefront analysis device with a pyramid or elongated Laser star), and
//  various systems configurations such as multi-conjugate AO.
//
//  COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
//  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//  See the GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License along with COMPASS.
//  If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>.
// -----------------------------------------------------------------------------

//! \file      sutra_groot.h
//! \ingroup   libsutra
//! \class     sutra_groot
//! \brief     this class provides the groot features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.4.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#ifndef _SUTRA_GROOT_H_
#define _SUTRA_GROOT_H_

#include <carma.h>
#include <carma_cublas.h>
#include <carma_host_obj.h>
#include <carma_obj.h>
#include <carma_sparse_obj.h>

class sutra_groot {
 public:
  carma_context *current_context;
  int device;

  int nactus;     // number of actuators
  int nssp;       // number of valid ssp
  int nlayers;    // number of atmos layers
  int npts;       // number of samples for aliasig computation
  float gsangle;  // Guide star angle [rad]
  float fc;       // DM cut-off frequency [m]
  float scale;
  float d;

  carma_host_obj<float> *h_vdt;      // v*dt/g [m/s]
  carma_host_obj<float> *h_Htheta;   // H*theta (theta = GS radial distance) [m]
  carma_host_obj<float> *h_L0;       // L0 [m]
  carma_host_obj<float> *h_winddir;  // wind directions [rad]
  carma_host_obj<float> *h_scale;    // r0**(-5/3) * frac * (lambda/2pi)**2
  carma_host_obj<float> *h_weights;  // Simpson method weights for integration
  // Covariance matrix estimation
  carma_obj<float> *d_Cerr;  // Residual error covariance matrix on DM actuators
  carma_obj<float>
      *d_TT;  // TT component of the residual error covariance matrix
  carma_obj<float> *d_CaXX;  // XX aliasing covariance matrix
  carma_obj<float> *d_CaYY;  // YY aliasing covariance matrix

  carma_obj<float>
      *d_TTPfilter;  // Tip-tilt and piston filter matrix (= Btt.dot(P))
  carma_obj<float> *d_pzt2tt;  // pzt to TT matrix
  carma_obj<float> *d_Nact;    // Coupling matrix
  carma_obj<float> *d_xpos;    // X-positions of DM actuators or ssp [m]
  carma_obj<float> *d_ypos;    // Y-positions of DM actuators or ssp [m]

  // Dphi lowpass
  carma_obj<float> *d_tab_int_x;  // Tabulated integral
  carma_obj<float> *d_tab_int_y;

 public:
  sutra_groot(carma_context *context, int device, int nactus, int nlayers,
              float gsangle, float *vdt, float *Htheta, float *L0,
              float *winddir, float *scale, float *pzt2tt, float *TTPfilter,
              float *Nact, float *xpos, float *ypos, float fc);

  sutra_groot(carma_context *context, int device, int nssp, float *weights,
              float scale, float *xpos, float *ypos, float fc, float d,
              int npts);
  ~sutra_groot();

  cublasHandle_t cublas_handle() { return current_context->get_cublasHandle(); }
  void init_common(carma_context *context, int device, float *xpos, float *ypos,
                   int N, float fc);
  int compute_Cerr();
  int compute_Calias();
};
template <class T_data>
int compute_Cerr_layer(T_data *Cerr, int N, T_data *tab_int_x,
                       T_data *tab_int_y, T_data *xpos, T_data *ypos,
                       T_data vdt, T_data Htheta, T_data L0, T_data fc,
                       T_data winddir, T_data gsangle, T_data scale, int Ntab,
                       carma_device *device);
template <class T_data>
int compute_Ca(T_data *CaXX, T_data *CaYY, int nssp, T_data *tab_int_x,
               T_data *tab_int_y, T_data *xpos, T_data *ypos, T_data offset,
               T_data d, T_data fc, T_data scale, T_data weight, int Ntab,
               carma_device *device);
template <class T_data>
int tab_u831J0(T_data *tab_int_x, T_data *tab_int_y, int N,
               carma_device *device);
template <class T_data>
void cumsum(T_data *odata, T_data *idata, int N);
template <class T_data>
int add_transpose(T_data *Cerr, int N, carma_device *device);

#endif
