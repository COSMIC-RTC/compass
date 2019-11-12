// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2019 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.
//  Distributed under GNU - LGPL
//
//  COMPASS is free software: you can redistribute it and/or modify it under the
//  terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or any later version.
//
//  COMPASS: End-to-end AO simulation tool using GPU acceleration
//  The COMPASS platform was designed to meet the need of high-performance for
//  the simulation of AO systems.
//
//  The final product includes a software package for simulating all the
//  critical subcomponents of AO, particularly in the context of the ELT and a
//  real-time core based on several control approaches, with performances
//  consistent with its integration into an instrument. Taking advantage of the
//  specific hardware architecture of the GPU, the COMPASS tool allows to
//  achieve adequate execution speeds to conduct large simulation campaigns
//  called to the ELT.
//
//  The COMPASS platform can be used to carry a wide variety of simulations to
//  both testspecific components of AO of the E-ELT (such as wavefront analysis
//  device with a pyramid or elongated Laser star), and various systems
//  configurations such as multi-conjugate AO.
//
//  COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY
//  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
//  FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
//  details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with COMPASS. If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>.
// -----------------------------------------------------------------------------

//! \file      sutra_centroider_maskedPix.h
//! \ingroup   libsutra
//! \class     sutra_centroider_maskedPix
//! \brief     this class provides the centroider_maskedPix features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.3.2
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#ifndef _SUTRA_CENTROIDER_MASKEDPIX_H_
#define _SUTRA_CENTROIDER_MASKEDPIX_H_

#include <sutra_centroider.h>
#include <sutra_wfs_pyr_pyrhr.h>
#include <string>

template <class Tin, class T>
class sutra_centroider_maskedPix : public sutra_centroider<Tin, T> {
 public:
  carma_obj<T> *d_selected_pix;
  carma_obj<T> *d_mask;

 public:
  sutra_centroider_maskedPix(carma_context *context, sutra_wfs *wfs,
                             long nvalid, long npupils, float offset,
                             float scale, bool filter_TT, int device);

  ~sutra_centroider_maskedPix();

  string get_type();

  int get_maskedPix(float *img, float *intensities, T *centroids, int *subindx,
                    int *subindy, int nvalid, int ns);
  int get_cog(float *img, float *intensities, T *centroids, int nvalid,
              int npix, int ntot);
  int get_cog(float *intensities, T *slopes, bool noise);
  int get_cog();
  int fill_selected_pix(carma_obj<T> *pix);
  int fill_mask();
};

void fill_intensities(float *intensities, float *img, int *subindx,
                      int *subindy, int ns, int nslopes, carma_device *device);
template <class T>
void getMaskedPix(T *centroids, T *ref, float *img, int *subindx, int *subindy,
                  float *psum, int ns, int nslopes, carma_device *device);
template <class T>
void pyr_fill_selected_pix(T *img, int img_sizex, T *pix, int *subindx, int *subindy,
                      int nvalid, carma_device *device);
template <class T>
void pyr_fill_mask(T *mask, int img_sizex, int *subindx, int *subindy,
                      int nvalid, carma_device *device);

#endif  // _SUTRA_CENTROIDER_MASKEDPIX_H_
