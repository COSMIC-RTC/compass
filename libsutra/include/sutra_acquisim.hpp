// This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
//
// COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU
// Lesser General Public License as published by the Free Software Foundation, either version 3 of
// the License, or any later version.
//
// COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
// even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License along with COMPASS. If
// not, see <https://www.gnu.org/licenses/>
//
//  Copyright (C) 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>

//! \file      sutra_acquisim.hpp
//! \ingroup   libsutra
//! \class     SutraAcquisim
//! \brief     this class provides the acquisition simulator to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#ifndef SUTRA_ACQUISIM_H_
#define SUTRA_ACQUISIM_H_

#include <sutra_sensors.hpp>
#include <sutra_wfs_sh.hpp>

class SutraAcquisim {
 public:
  int32_t device;
  string type;
  int64_t nxsub;
  int64_t nvalid;
  int64_t npix;

  CarmaContext *current_context;

  SutraWfsSH *wfs;
  CarmaObj<int32_t> *d_validsubsx;
  CarmaObj<int32_t> *d_validsubsy;

 public:
  SutraAcquisim(SutraSensors *sensors, int32_t wfs_num);
  SutraAcquisim(const SutraAcquisim &acquisim);
  ~SutraAcquisim();

  int32_t set_validsubs(int64_t nvalid, int32_t *validsubsx, int32_t *validsubsy);

  int32_t comp_image_tele(int64_t *dims, float *bimage);
  int32_t comp_image(int64_t *dims, float *bimage);
  int32_t comp_image(int64_t *dims, float *bimage, CarmaObj<float> *d_bincube);
  int32_t comp_image_2D(int64_t *dims, float *bimage, int32_t *num_ssp);

 private:
};

// General utilities
template <class T>
int32_t fillbincube_2D(T *bimage, T *bcube, int32_t npix, int32_t nsub, int32_t *valid);

template <class T>
int32_t fillbincube(T *bimage, T *bcube, int32_t npix, int32_t nsub, int32_t Nsub, int32_t *ivalid,
                int32_t *jvalid, CarmaDevice *device);

template <class T>
int32_t fillbincube_async(CarmaStreams *streams, CarmaObj<T> *bimage,
                      CarmaObj<T> *bcube, int32_t npix, int32_t nsub, int32_t Nsub,
                      int32_t *ivalid, int32_t *jvalid, CarmaDevice *device);

template <class T>
int32_t fillbincube_async(CarmaHostObj<T> *image_telemetry, T *bimage, T *bcube,
                      int32_t npix, int32_t nsub, int32_t Nsub, int32_t *ivalid, int32_t *jvalid,
                      int32_t nim, CarmaDevice *device);

#endif /* SUTRA_ACQUISIM_H_ */
