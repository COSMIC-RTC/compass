// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_centroider_corr.h
//! \ingroup   libsutra
//! \class     SutraCentroiderCorr
//! \brief     this class provides the centroider_corr features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24


#ifndef _SUTRA_CENTROIDER_CORR_H_
#define _SUTRA_CENTROIDER_CORR_H_

#include <sutra_centroider.h>

template <class Tin, class T>
class SutraCentroiderCorr : public SutraCentroider<Tin, T> {
 public:
  int interp_sizex;
  int interp_sizey;
  CarmaObj<cuFloatComplex> *d_corrfnct;
  CarmaObj<cuFloatComplex> *d_corrspot;
  CarmaObj<T> *d_corrnorm;
  CarmaObj<int> *d_corrmax;
  CarmaObj<T> *d_corr;
  CarmaObj<T> *d_interpmat;

 public:
  SutraCentroiderCorr(CarmaContext *context, SutraWfs *wfs, long nvalid,
                        float offset, float scale, bool filter_TT, int device);
  SutraCentroiderCorr(const SutraCentroiderCorr &centroider);
  ~SutraCentroiderCorr();

  string get_type();
  int fill_bincube(T *img);

  int init_corr(int isizex, int isizey, T *interpmat);
  int load_corr(T *corr, T *corr_norm, int ndim);

  int get_cog(float *cube, float *intensities, T *centroids, int nvalid,
              int npix, int ntot, cudaStream_t stream=0);
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
             CarmaDevice *device);

template <class T>
int correl(T *d_odata, T *d_idata, int N, CarmaDevice *device);

template <class Tcu, class T>
int roll2real(T *d_odata, Tcu *d_idata, int n, int Npix, int N,
              CarmaDevice *device);

template <class T>
int corr_norm(T *d_odata, T *d_idata, int Npix, int N, CarmaDevice *device);

#endif  // _SUTRA_CENTROIDER_CORR_H_
