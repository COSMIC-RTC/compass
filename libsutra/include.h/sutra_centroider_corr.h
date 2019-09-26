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

//! \file      sutra_centroider_corr.h
//! \ingroup   libsutra
//! \class     sutra_centroider_corr
//! \brief     this class provides the centroider_corr features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.3.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License


#ifndef _SUTRA_CENTROIDER_CORR_H_
#define _SUTRA_CENTROIDER_CORR_H_

#include <sutra_centroider.h>

template <class Tin, class T>
class sutra_centroider_corr : public sutra_centroider<Tin, T> {
 public:
  int npix;
  int interp_sizex;
  int interp_sizey;
  carma_obj<cuFloatComplex> *d_corrfnct;
  carma_obj<cuFloatComplex> *d_corrspot;
  carma_obj<T> *d_corrnorm;
  carma_obj<int> *d_corrmax;
  carma_obj<T> *d_corr;
  carma_obj<T> *d_interpmat;

 public:
  sutra_centroider_corr(carma_context *context, sutra_wfs *wfs, long nvalid,
                        float offset, float scale, bool filter_TT, int device);
  sutra_centroider_corr(const sutra_centroider_corr &centroider);
  ~sutra_centroider_corr();

  string get_type();
  int fill_bincube(T *img);

  int init_corr(int isizex, int isizey, T *interpmat);
  int load_corr(T *corr, T *corr_norm, int ndim);

  int set_npix(int npix);

  int get_cog(float *cube, float *intensities, T *centroids, int nvalid,
              int npix, int ntot);
  int get_cog(float *intensities, T *slopes, bool noise);
  int get_cog();
};

template <class T>
void subap_sortmaxi(int threads, int blocks, T *d_idata, int *values, int nmax,
                    int offx, int offy, int npix, int Npix);
template <class T>
void subap_pinterp(int threads, int blocks, T *d_idata, int *values,
                   T *d_centroids, T *d_matinterp, int sizex, int sizey,
                   int nvalid, int Npix, float scale, float offset);

template <class Tcu, class T>
int fillcorr(Tcu *d_out, T *d_in, int npix_in, int npix_out, int N, int nvalid,
             carma_device *device);

template <class T>
int correl(T *d_odata, T *d_idata, int N, carma_device *device);

template <class Tcu, class T>
int roll2real(T *d_odata, Tcu *d_idata, int n, int Npix, int N,
              carma_device *device);

template <class T>
int corr_norm(T *d_odata, T *d_idata, int Npix, int N, carma_device *device);

#endif  // _SUTRA_CENTROIDER_CORR_H_
