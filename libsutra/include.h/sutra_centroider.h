// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_centroider.h
//! \ingroup   libsutra
//! \class     SutraCentroider
//! \brief     this class provides the centroider features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.3.0
//! \date      2022/01/24
//! \copyright GNU Lesser General Public License


#ifndef _SUTRA_CENTROIDER_H_
#define _SUTRA_CENTROIDER_H_

#include <sutra_acquisim.h>
#include <sutra_wfs.h>
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
  int device;
  SutraWfs *wfs;
  int nvalid;
  int nslopes;
  int npix;
  int nxsub;
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
  CarmaObj<int> *d_lutPix;
  CarmaObj<int> *d_validx;
  CarmaObj<int> *d_validy;
  CarmaObj<int> *d_validMask;

  CarmaObj<float> *d_centro_filtered;
  CarmaObj<float> *d_ref_Tip;
  CarmaObj<float> *d_ref_Tilt;
  CarmaObj<float> *d_TT_slopes;

protected:
  SutraCentroider(CarmaContext *context, SutraWfs *wfs, long nvalid,
                   float offset, float scale, bool filter_TT, int device);

private:
  template <typename Q = Tout>
  typename std::enable_if<std::is_same<Q, float>::value, int>::type
  apply_TT_filter_impl(Tout *centroids, std::true_type);
  int apply_TT_filter_impl(Tout *centroids, std::false_type);

public:
  virtual ~SutraCentroider();
  int set_scale(float scale);
  int set_offset(float offset);
  int set_dark(float *dark, int n);
  int set_flat(float *flat, int n);
  int set_lutPix(int *lutPix, int n);
  int init_calib(int n, int m);
  int init_roi(int N);
  int set_centroids_ref(float *centroids_ref);
  int calibrate_img_validPix() { return calibrate_img_validPix(0); };
  int calibrate_img_validPix(cudaStream_t stream);
  int calibrate_img() { return calibrate_img(0); };
  int calibrate_img(cudaStream_t stream);
  int load_validpos(int *ivalid, int *jvalid, int N);
  int set_npix(int npix);
  int set_nxsub(int nxsub);
  int load_img(Tin *img, int n);
  int load_img(Tin *img, int n, int location);
  int load_img(Tin *img, int m, int n, int location);
  int load_img(CarmaObj<Tin> *img);
  int init_img_raw(int m, int n);
  int get_validMask();
  bool is_type(string typec) { return (typec.compare(get_type()) == 0); }
  int init_TT_filter();
  int apply_TT_filter(Tout *centroids);

  virtual string get_type() = 0;

  virtual int get_cog(float *img, float *intensities, Tout *centroids,
                      int nvalid, int npix, int ntot,
                      cudaStream_t stream = 0) = 0;
  virtual int get_cog(float *intensities, Tout *slopes, bool noise) = 0;
  virtual int get_cog() = 0;
};

template <class Tin>
int calibration_validPix_sh(int npix, int size, int blocks, Tin *img_raw,
                            float *img_cal, float *dark, float *flat,
                            int *lutPix, int *validx, int *validy,
                            CarmaDevice *device, cudaStream_t stream);

template <class Tin>
int calibration_validPix_pyr(Tin *img_raw, float *img_cal, float *dark,
                             float *flat, int *lutPix, int *validx, int *validy,
                             int nvalid, int img_sizex, CarmaDevice *device,
                             cudaStream_t stream = 0);

template <class Tin>
int calibration(Tin *img_raw, float *img_cal, float *dark, float *flat,
                int *lutPix, int N, CarmaDevice *device,
                cudaStream_t stream = 0);

template <typename T>
int convert_centro(T *d_odata, T *d_idata, float offset, float scale, int N,
                   CarmaDevice *device);
int fill_validMask(int size, int npix, int blocks, int *d_validMask,
                   int *validx, int *validy, CarmaDevice *device);

#endif  // _SUTRA_CENTROIDER_H_
