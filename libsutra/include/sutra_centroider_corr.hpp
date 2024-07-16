// This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
//
// COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
// General Public License as published by the Free Software Foundation, either version 3 of the 
// License, or any later version.
//
// COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
// See the GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License along with COMPASS. 
// If not, see <https://www.gnu.org/licenses/>
//
//  Copyright (C) 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>

//! \file      sutra_centroider_corr.hpp
//! \ingroup   libsutra
//! \class     SutraCentroiderCorr
//! \brief     this class provides the centroider_corr features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24


#ifndef _SUTRA_CENTROIDER_CORR_H_
#define _SUTRA_CENTROIDER_CORR_H_

#include <sutra_centroider.hpp>

template <class Tin, class Tcomp>
class SutraCentroiderCorr : public SutraCentroider<Tin, Tcomp> {
 public:
  int32_t interp_sizex;
  int32_t interp_sizey;
  CarmaObj<cuFloatComplex> *d_corrfnct;
  CarmaObj<cuFloatComplex> *d_corrspot;
  CarmaObj<Tcomp> *d_corrnorm;
  CarmaObj<int32_t> *d_corrmax;
  CarmaObj<Tcomp> *d_corr;
  CarmaObj<Tcomp> *d_interpmat;

 public:
  SutraCentroiderCorr(CarmaContext *context, SutraWfs *wfs, int64_t nvalid,
                        float offset, float scale, bool filter_TT, int32_t device);
  SutraCentroiderCorr(const SutraCentroiderCorr &centroider);
  ~SutraCentroiderCorr();

  string get_type();
  int32_t fill_bincube(Tcomp *img);

  int32_t init_corr(int32_t isizex, int32_t isizey, Tcomp *interpmat);
  int32_t load_corr(Tcomp *corr, Tcomp *corr_norm, int32_t ndim);

  int32_t get_cog(float *cube, float *intensities, Tcomp *centroids, int32_t nvalid,
              int32_t npix, int32_t ntot, cudaStream_t stream=0);
  int32_t get_cog(float *intensities, Tcomp *slopes, bool noise);
  int32_t get_cog();
};

template <class Tcomp>
void subap_sortmaxi(int32_t threads, int32_t blocks, Tcomp *d_idata, int32_t *values, int32_t nmax,
                    int32_t offx, int32_t offy, int32_t npix, int32_t Npix);
template <class Tcomp>
void subap_pinterp(int32_t threads, int32_t blocks, Tcomp *d_idata, int32_t *values,
                   Tcomp *d_centroids, Tcomp *d_matinterp, int32_t sizex, int32_t sizey,
                   int32_t nvalid, int32_t Npix, float scale, float offset);

template <class Tcu, class Tcomp>
int32_t fillcorr(Tcu *d_out, Tcomp *d_in, int32_t npix_in, int32_t npix_out, int32_t N, int32_t nvalid,
             CarmaDevice *device);

template <class Tcomp>
int32_t correl(Tcomp *d_odata, Tcomp *d_idata, int32_t N, CarmaDevice *device);

template <class Tcu, class Tcomp>
int32_t roll2real(Tcomp *d_odata, Tcu *d_idata, int32_t n, int32_t Npix, int32_t N,
              CarmaDevice *device);

template <class Tcomp>
int32_t corr_norm(Tcomp *d_odata, Tcomp *d_idata, int32_t Npix, int32_t N, CarmaDevice *device);

#endif  // _SUTRA_CENTROIDER_CORR_H_
