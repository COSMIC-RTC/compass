// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

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
