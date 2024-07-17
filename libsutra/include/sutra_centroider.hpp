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

//! \file      sutra_centroider.hpp
//! \ingroup   libsutra
//! \class     SutraCentroider
//! \brief     this class provides the centroider features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24


#ifndef _SUTRA_CENTROIDER_H_
#define _SUTRA_CENTROIDER_H_

#include <sutra_acquisim.hpp>
#include <sutra_wfs.hpp>
#include <string>

enum class SlopeOrder {
    untied     = 0, // x,x,x,y,y,y
    interlaced = 1  // x,y,x,y,x,y
};

inline SlopeOrder slope_order(std::size_t value) {
    return static_cast<SlopeOrder>(value);
}
template <class Tin, class Tout>
class SutraCentroider {
public:
  int32_t device;
  SutraWfs *wfs;
  int32_t nvalid;
  int32_t nslopes;
  int32_t npix;
  int32_t nxsub;
  bool filter_TT;

  float offset;
  float scale;

  SlopeOrder slope_order = SlopeOrder::untied;

  CarmaContext *current_context;

  CarmaObj<Tout> *d_bincube;
  CarmaObj<float> *d_intensities;
  CarmaObj<Tout> *d_centroids_ref;  // ref centroids
  CarmaObj<float> *d_img;
  CarmaObj<Tin> *d_img_raw;
  CarmaObj<float> *d_dark;
  CarmaObj<float> *d_flat;
  CarmaObj<int32_t> *d_lutPix;
  CarmaObj<int32_t> *d_validx;
  CarmaObj<int32_t> *d_validy;
  CarmaObj<int32_t> *d_validMask;

  CarmaObj<float> *d_centro_filtered;
  CarmaObj<float> *d_ref_Tip;
  CarmaObj<float> *d_ref_Tilt;
  CarmaObj<float> *d_TT_slopes;

protected:
  SutraCentroider(CarmaContext *context, SutraWfs *wfs, int64_t nvalid,
                   float offset, float scale, bool filter_TT, int32_t device);

public:
  virtual ~SutraCentroider();
  int32_t set_scale(float scale);
  int32_t set_offset(float offset);
  int32_t set_dark(float *dark, int32_t n);
  int32_t set_flat(float *flat, int32_t n);
  int32_t set_lutPix(int32_t *lutPix, int32_t n);
  int32_t init_calib(int32_t n, int32_t m);
  int32_t init_roi(int32_t N);
  int32_t set_centroids_ref(float *centroids_ref);
  int32_t calibrate_img_validPix() { return calibrate_img_validPix(0); };
  int32_t calibrate_img_validPix(cudaStream_t stream);
  int32_t calibrate_img() { return calibrate_img(0); };
  int32_t calibrate_img(cudaStream_t stream);
  int32_t load_validpos(int32_t *ivalid, int32_t *jvalid, int32_t N);
  int32_t set_npix(int32_t npix);
  int32_t set_nxsub(int32_t nxsub);
  int32_t load_img(Tin *img, int32_t n);
  int32_t load_img(Tin *img, int32_t n, int32_t location);
  int32_t load_img(Tin *img, int32_t m, int32_t n, int32_t location);
  int32_t load_img(CarmaObj<Tin> *img);
  int32_t init_img_raw(int32_t m, int32_t n);
  int32_t get_validMask();
  bool is_type(string typec) { return (typec.compare(get_type()) == 0); }
  int32_t init_TT_filter();
  int32_t apply_TT_filter(Tout *centroids);

  virtual string get_type() = 0;

  virtual int32_t get_cog(float *img, float *intensities, Tout *centroids,
                      int32_t nvalid, int32_t npix, int32_t ntot,
                      cudaStream_t stream = 0) = 0;
  virtual int32_t get_cog(float *intensities, Tout *slopes, bool noise) = 0;
  virtual int32_t get_cog() = 0;
};

template <class Tin>
int32_t calibration_validPix_sh(int32_t npix, int32_t size, int32_t blocks, Tin *img_raw,
                            float *img_cal, float *dark, float *flat,
                            int32_t *lutPix, int32_t *validx, int32_t *validy,
                            CarmaDevice *device, cudaStream_t stream);

template <class Tin>
int32_t calibration_validPix_pyr(Tin *img_raw, float *img_cal, float *dark,
                             float *flat, int32_t *lutPix, int32_t *validx, int32_t *validy,
                             int32_t nvalid, int32_t img_sizex, CarmaDevice *device,
                             cudaStream_t stream = 0);

template <class Tin>
int32_t calibration(Tin *img_raw, float *img_cal, float *dark, float *flat,
                int32_t *lutPix, int32_t N, CarmaDevice *device,
                cudaStream_t stream = 0);

template <typename T>
int32_t convert_centro(T *d_odata, T *d_idata, float offset, float scale, int32_t N,
                   CarmaDevice *device);
int32_t fill_validMask(int32_t size, int32_t npix, int32_t blocks, int32_t *d_validMask,
                   int32_t *validx, int32_t *validy, CarmaDevice *device);

#endif  // _SUTRA_CENTROIDER_H_
