// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_acquisim.hpp
//! \ingroup   libsutra
//! \class     SutraAcquisim
//! \brief     this class provides the acquisition simulator to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
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
