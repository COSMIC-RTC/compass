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

//! \file      sutra_centroider.h
//! \ingroup   libsutra
//! \class     SutraCentroider
//! \brief     this class provides the centroider features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.0.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License


#ifndef _SUTRA_CENTROIDER_H_
#define _SUTRA_CENTROIDER_H_

#include <sutra_acquisim.h>
#include <sutra_wfs.h>
#include <string>

template <class Tin, class Tout>
class SutraCentroider {
public:
  int device;
  SutraWfs *wfs;
  int nvalid;
  int nslopes;
  int npix;
  int nxsub;
  bool filter_TT;

  float offset;
  float scale;

  CarmaContext *current_context;

  CarmaObj<Tout> *d_bincube;
  CarmaObj<float> *d_intensities;
  CarmaObj<Tout> *d_centroids_ref;  // ref centroids
  CarmaObj<float> *d_img;
  CarmaObj<Tin> *d_img_raw;
  CarmaObj<float> *d_dark;
  CarmaObj<float> *d_flat;
  CarmaObj<int> *d_lutPix;
  CarmaObj<int> *d_validx;
  CarmaObj<int> *d_validy;
  CarmaObj<int> *d_validMask;

  CarmaObj<float> *d_centro_filtered;
  CarmaObj<float> *d_ref_Tip;
  CarmaObj<float> *d_ref_Tilt;
  CarmaObj<float> *d_TT_slopes;

protected:
  SutraCentroider(CarmaContext *context, SutraWfs *wfs, long nvalid,
                   float offset, float scale, bool filter_TT, int device);

private:
  template <typename Q = Tout>
  typename std::enable_if<std::is_same<Q, float>::value, int>::type
  apply_TT_filter_impl(Tout *centroids, std::true_type);
  int apply_TT_filter_impl(Tout *centroids, std::false_type);

public:
  virtual ~SutraCentroider();
  int set_scale(float scale);
  int set_offset(float offset);
  int set_dark(float *dark, int n);
  int set_flat(float *flat, int n);
  int set_lutPix(int *lutPix, int n);
  int init_calib(int n, int m);
  int init_roi(int N);
  int set_centroids_ref(float *centroids_ref);
  int calibrate_img() {return calibrate_img(0);};
  int calibrate_img(cudaStream_t stream);
  int load_validpos(int *ivalid, int *jvalid, int N);
  int set_npix(int npix);
  int set_nxsub(int nxsub);
  int load_img(Tin *img, int n);
  int load_img(Tin *img, int n, int location);
  int load_img(Tin *img, int m, int n, int location);
  int load_img(CarmaObj<Tin> *img);
  int init_img_raw(int m, int n);
  int get_validMask();
  bool is_type(string typec) { return (typec.compare(get_type()) == 0); }
  int init_TT_filter();
  int apply_TT_filter(Tout *centroids);

  virtual string get_type() = 0;

  virtual int get_cog(float *img, float *intensities, Tout *centroids,
                      int nvalid, int npix, int ntot, cudaStream_t stream=0) = 0;
  virtual int get_cog(float *intensities, Tout *slopes, bool noise) = 0;
  virtual int get_cog() = 0;
};
template <class Tin>
int calibration(Tin *img_raw, float *img_cal, float *dark, float *flat, int *lutPix, int N,
                CarmaDevice *device, cudaStream_t stream=0);

template <typename T>
int convert_centro(T *d_odata, T *d_idata, float offset, float scale, int N,
                   CarmaDevice *device);
int fill_validMask(int size, int npix, int blocks, int *d_validMask,
                   int *validx, int *validy, CarmaDevice *device);

#endif  // _SUTRA_CENTROIDER_H_
