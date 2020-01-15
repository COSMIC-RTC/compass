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

//! \file      sutra_centroider_wcog.h
//! \ingroup   libsutra
//! \class     sutra_centroider_wcog
//! \brief     this class provides the centroider_wcog features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.4.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License


#ifndef _SUTRA_CENTROIDER_WCOG_H_
#define _SUTRA_CENTROIDER_WCOG_H_

#include <sutra_centroider.h>

template <class Tin, class T>
class sutra_centroider_wcog : public sutra_centroider<Tin, T> {
 public:
  int npix;
  carma_obj<float> *d_weights;

 public:
  sutra_centroider_wcog(carma_context *context, sutra_wfs *wfs, long nvalid,
                        float offset, float scale, bool filter_TT, int device);
  sutra_centroider_wcog(const sutra_centroider_wcog &centroider);
  ~sutra_centroider_wcog();

  string get_type();

  int set_npix(int npix);
  int init_weights();
  int load_weights(float *weights, int ndim);

  int get_cog(float *cube, float *intensities, T *centroids, int nvalid,
              int npix, int ntot);
  int get_cog(float *intensities, T *slopes, bool noise);
  int get_cog();
};

template <class T>
void get_centroids(int size, int threads, int blocks, int n, float *d_idata,
                   T *d_odata, T *ref, int *validx, int *validy,
                   float *intensities, float *d_weights, float scale,
                   float offset, carma_device *device);

template <class T>
int fillweights(T *d_out, T *d_in, int npix, int N, carma_device *device);
#endif  // _SUTRA_CENTROIDER_WCOG_H_
