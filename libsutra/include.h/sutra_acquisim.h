// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
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

//! \file      sutra_acquisim.h
//! \ingroup   libsutra
//! \class     SutraAcquisim
//! \brief     this class provides the acquisition simulator to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.3.0
//! \date      2022/01/24
//! \copyright GNU Lesser General Public License

#ifndef SUTRA_ACQUISIM_H_
#define SUTRA_ACQUISIM_H_

#include <sutra_sensors.h>
#include <sutra_wfs_sh.h>

class SutraAcquisim {
 public:
  int device;
  string type;
  long nxsub;
  long nvalid;
  long npix;

  CarmaContext *current_context;

  SutraWfsSH *wfs;
  CarmaObj<int32_t> *d_validsubsx;
  CarmaObj<int32_t> *d_validsubsy;

 public:
  SutraAcquisim(SutraSensors *sensors, int wfs_num);
  SutraAcquisim(const SutraAcquisim &acquisim);
  ~SutraAcquisim();

  int set_validsubs(int64_t nvalid, int32_t *validsubsx, int32_t *validsubsy);

  int comp_image_tele(long *dims, float *bimage);
  int comp_image(long *dims, float *bimage);
  int comp_image(long *dims, float *bimage, CarmaObj<float> *d_bincube);
  int comp_image_2D(long *dims, float *bimage, int *num_ssp);

 private:
};

// General utilities
template <class T>
int fillbincube_2D(T *bimage, T *bcube, int npix, int nsub, int *valid);

template <class T>
int fillbincube(T *bimage, T *bcube, int npix, int nsub, int Nsub, int *ivalid,
                int *jvalid, CarmaDevice *device);

template <class T>
int fillbincube_async(CarmaStreams *streams, CarmaObj<T> *bimage,
                      CarmaObj<T> *bcube, int npix, int nsub, int Nsub,
                      int *ivalid, int *jvalid, CarmaDevice *device);

template <class T>
int fillbincube_async(CarmaHostObj<T> *image_telemetry, T *bimage, T *bcube,
                      int npix, int nsub, int Nsub, int *ivalid, int *jvalid,
                      int nim, CarmaDevice *device);

#endif /* SUTRA_ACQUISIM_H_ */
