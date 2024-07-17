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

//! \file      sutra_wfs.hpp
//! \ingroup   libsutra
//! \class     SutraWfs
//! \brief     this class provides the wfs features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#ifndef _SUTRA_WFS_H_
#define _SUTRA_WFS_H_

#include <carma_utils.hpp>
#include <sutra_lgs.hpp>
#include <sutra_phase.hpp>
#include <sutra_target.hpp>
#include <sutra_telemetry.hpp>
#include <sutra_telescope.hpp>
#include <sutra_utils.hpp>

#include <map>
#include <vector>
//#include <sutra_slopes.hpp>

using std::map;
using std::string;

class SutraWfs {
 public:
  int32_t device;
  string type;
  int64_t nxsub;
  int64_t nvalid;
  int64_t npix;
  int64_t nrebin;
  int64_t nfft;
  int64_t ntot;
  int64_t npup;
  int64_t nphase;
  int64_t nmaxhr;
  int64_t nffthr;
  float subapd;
  float nphot;
  float nphot4imat;
  float noise;
  bool lgs;
  bool kernconv;
  bool roket;
  bool is_low_order;

  bool fakecam;
  int32_t max_flux_per_pix;
  int32_t max_pix_value;

  cufftHandle *campli_plan;
  cufftHandle *fttotim_plan;
  CarmaObj<cuFloatComplex> *d_ftkernel;
  CarmaObj<cuFloatComplex> *d_camplipup;
  CarmaObj<cuFloatComplex> *d_camplifoc;
  CarmaObj<cuFloatComplex> *d_fttotim;

  CarmaObj<float> *d_pupil;
  CarmaObj<float> *d_bincube;
  CarmaObj<float> *d_binimg;
  CarmaObj<float> *d_binimg_notnoisy;
  CarmaObj<float> *d_intensities;
  CarmaObj<float> *d_offsets;
  CarmaObj<float> *d_fluxPerSub;
  CarmaObj<float> *d_sincar;
  CarmaObj<int32_t> *d_hrmap;
  CarmaObj<uint16_t> *d_camimg;
  CarmaObj<float> *d_dark;
  CarmaObj<float> *d_flat;

  CarmaObj<float> *d_slopes;

  CarmaHostObj<float> *image_telemetry;

  SutraSource *d_gs;
  std::vector<CarmaObj<float> *> d_pupil_ngpu;

  CarmaStreams *streams;
  int32_t nstreams;

  CarmaObj<int32_t> *d_phasemap;
  CarmaObj<float> *d_ttprojmat;
  CarmaObj<float> *d_ttprojvec;
  CarmaObj<int32_t> *d_validsubsx;  // nvalid
  CarmaObj<int32_t> *d_validsubsy;  // nvalid
  CarmaObj<float> *d_submask; //field stop

  CarmaContext *current_context;

  /// MPI stuff
  int32_t offset;
  int32_t nvalid_tot;
  int32_t rank;
  int32_t *displ_bincube;
  int32_t *count_bincube;

 public:
  virtual ~SutraWfs(){};

  int32_t wfs_initgs(CarmaObj<float> *d_lgskern,
                 CarmaObj<cuFloatComplex> *d_ftlgskern,
                 map<vector<int32_t>, cufftHandle *> ftlgskern_plans, float xpos,
                 float ypos, float lambda, float mag, float zerop, int64_t size,
                 float noise, int64_t seed, float G, float thetaML, float dx,
                 float dy);
  int32_t set_pupil(float *pupil);
  int32_t set_binimg(float *binimg, int32_t nElem);
  int32_t set_dark(float *dark, int32_t nElem);
  int32_t set_flat(float *flat, int32_t nElem);
  int32_t set_fakecam(bool fakecam);
  int32_t set_max_flux_per_pix(int32_t max_flux_per_pix);
  int32_t set_max_pix_value(int32_t max_pix_value);

  int32_t load_kernels(float *lgskern);
  int32_t sensor_trace(SutraAtmos *yatmos);
  int32_t sensor_trace(SutraDms *ydm, int32_t rst);
  int32_t sensor_trace(SutraAtmos *atmos, SutraDms *ydms);
  int32_t sensor_trace(int32_t rst);
  int32_t slopes_geom(float *slopes, int32_t type = 0);
  int32_t slopes_geom(int32_t type = 0);

  virtual int32_t fill_binimage(int32_t async) = 0;
  virtual int32_t comp_image(bool noise = true) = 0;

  virtual int32_t define_mpi_rank(int32_t rank, int32_t size) = 0;
  virtual int32_t allocate_buffers(
      map<vector<int32_t>, cufftHandle *> campli_plans,
      map<vector<int32_t>, cufftHandle *> fttotim_plans) = 0;
  int32_t set_noise(float noise, int64_t seed);

 protected:
  virtual int32_t comp_generic() = 0;
  SutraWfs(CarmaContext *context, SutraTelescope *d_tel,
            CarmaObj<cuFloatComplex> *d_camplipup,
            CarmaObj<cuFloatComplex> *d_camplifoc,
            CarmaObj<cuFloatComplex> *d_fttotim, string type, int64_t nxsub,
            int64_t nvalid, int64_t npix, int64_t nphase, int64_t nrebin, int64_t nfft,
            int64_t ntot, int64_t npup, float pdiam, float nphotons, float nphot4imat,
            int32_t lgs, bool fakecam, int32_t max_flux_per_pix, int32_t max_pix_value,
            bool is_low_order, bool roket, int32_t device);
};

// General utilities
int32_t fillcamplipup(cuFloatComplex *amplipup, float *phase, float *offset,
                  float *mask, float scale, int32_t *istart, int32_t *jstart,
                  int32_t *ivalid, int32_t *jvalid, int32_t nphase, int32_t npup, int32_t Nfft,
                  int32_t Ntot, CarmaDevice *device, int32_t offset_phase);
int32_t fillcamplipup(cuFloatComplex *amplipup, cuFloatComplex *phase, float *offset,
                  int32_t *istart, int32_t *jstart,
                  int32_t *ivalid, int32_t *jvalid, int32_t nphase, int32_t N, int32_t Nfft,
                  int32_t Ntot, CarmaDevice *device);
int32_t fillfsamplipup(cuFloatComplex *d_odata, float *idata, float *mask, float scale,
                    int32_t Nfft, int32_t N, CarmaDevice *device);
int32_t indexfill(cuFloatComplex *d_odata, cuFloatComplex *d_idata, int32_t *indx,
              int32_t nx, int32_t Nx, int32_t N, CarmaDevice *device);
int32_t fillbincube(float *bcube, cuFloatComplex *hrimage, int32_t *indxpix, int32_t Nfft,
                int32_t Npix, int32_t Nrebin, int32_t Nsub, CarmaDevice *device);
int32_t fillbincube_async(CarmaStreams *streams, float *bcube,
                      cuFloatComplex *hrimage, int32_t *indxpix, int32_t Nfft, int32_t Npix,
                      int32_t Nrebin, int32_t Nsub, CarmaDevice *device);
int32_t fillbinimg(float *bimage, float *bcube, int32_t npix, int32_t nsub, int32_t Nsub,
               int32_t *ivalid, int32_t *jvalid, bool add, CarmaDevice *device);
int32_t fillbinimg_async(CarmaStreams *streams, CarmaObj<float> *bimage,
                     CarmaObj<float> *bcube, int32_t npix, int32_t nsub, int32_t Nsub,
                     int32_t *ivalid, int32_t *jvalid, bool add, CarmaDevice *device);
int32_t fillbinimg_async(CarmaHostObj<float> *image_telemetry, float *bimage,
                     float *bcube, int32_t npix, int32_t nsub, int32_t Nsub, int32_t *ivalid,
                     int32_t *jvalid, int32_t nim, bool add, CarmaDevice *device);
int32_t convolve_cube(cuFloatComplex *d_odata, cuFloatComplex *d_idata, int32_t N,
                  int32_t n, CarmaDevice *device);
template <class T>
int32_t digitalize(T *camimg, float *binimg, float *dark, float *flat,
               int32_t max_flux_per_pix, int32_t max_pix_value, int32_t N, CarmaDevice *device);
// CUDA templates
// this is for cog
template <class T>
void subap_reduce(int32_t size, int32_t threads, int32_t blocks, T *d_idata, T *d_odata,
                  CarmaDevice *device);

template <class T>
void subap_reduce_async(int32_t threads, int32_t blocks, CarmaStreams *streams,
                        T *d_idata, T *d_odata);
// this is for tcog
template <class T>
void subap_reduce(int32_t size, int32_t threads, int32_t blocks, T *d_idata, T *d_odata,
                  T thresh, CarmaDevice *device);
// this is for tcog_new
template <class T>
void subap_reduce_new(int32_t size, int32_t threads, int32_t blocks, T *d_idata, T *d_odata,
                      T thresh, CarmaDevice *device);
// this is for wcog
template <class T>
void subap_reduce(int32_t size, int32_t threads, int32_t blocks, T *d_idata, T *d_odata,
                  T *weights, CarmaDevice *device);

template <class T>
void phase_reduce(int32_t threads, int32_t blocks, T *d_idata, T *d_odata, int32_t *indx,
                  T alpha);

template <class T>
void phase_derive(int32_t size, int32_t threads, int32_t blocks, int32_t n, T *d_idata,
                  T *d_odata, int32_t *indx, T *mask, T alpha, float *fluxPerSub);

template <class T>
void phase_project(int32_t nphase, int32_t nvalid, T *d_idata, T *d_odata, int32_t *indx,
                   T *d_ttprojmat, T *d_ttprojvec, CarmaDevice *device);

template <class Tout, class Tin>
void pyr_getpup(Tout *d_odata, Tin *d_idata, Tout *d_offsets, Tin *d_pup,
                int32_t np, float lambda, CarmaDevice *device);

template <class Tout, class Tin>
void pyr_getpup(Tout *d_odata, Tin *d_idata, Tin *d_pup, int32_t np, int32_t N,
                float lambda, float cx, float cy, CarmaDevice *device);

template <class T>
void pyr_rollmod(T *d_odata, T *d_idata, T *d_mask, float cx, float cy, int32_t np,
                 int32_t ns, CarmaDevice *device);

template <class T>
void pyr_fillbin(T *d_odata, T *d_idata, int32_t nrebin, int32_t np, int32_t ns, int32_t nim,
                 CarmaDevice *device);

template <class T>
int32_t pyr_fillbinimg(T *bimage, const T *bcube, const int32_t nxsub, const bool add,
                   CarmaDevice *device);

template <class T>
int32_t pyr_fillbinimg(T *oimage, const T *image, const int32_t n, const int32_t N,
                   const int32_t rebin, const bool add, CarmaDevice *device);

template <class Tin, class Tout>
void pyr_abs2(Tout *d_odata, Tin *d_idata, Tout fact, int32_t ns, int32_t nim,
              CarmaDevice *device);

template <class Tout, class Tin>
void apply_submask(Tout *d_odata, Tin *d_mask, int32_t n, CarmaDevice *device);

template <class T>
void pyr_submaskpyr(T *d_odata, T *d_mask, int32_t n, CarmaDevice *device);

template <class T>
void pyr_submaskpyr(T *d_odata, float *d_mask, int32_t n, CarmaDevice *device);

template <class Tout, class Tin>
void pyr_abs(Tout *d_odata, Tin *d_idata, int32_t ns, int32_t nim,
             CarmaDevice *device);

template <class Tout, class Tin>
void pyr_submask3d(Tout *d_odata, Tin *d_mask, int32_t n, int32_t nim,
                   CarmaDevice *device);

template <class T>
void pyr_intensities(T *d_odata, T *d_idata, int32_t *subindx, int32_t *subindy, int32_t ns,
                     int32_t nvalid, int32_t nim, CarmaDevice *device, cudaStream_t stream=0);

template <class T>
void pyr_intensities(T *d_odata, T *d_idata, int32_t *subindx, int32_t *subindy, int32_t ns,
                     int32_t nvalid, CarmaDevice *device, cudaStream_t stream=0);

template <class T>
void pyr_fact(T *d_data, T fact, int32_t n, int32_t nim, CarmaDevice *device);

void pyr_fact(cuFloatComplex *d_data, float fact, int32_t n, int32_t nim,
              CarmaDevice *device);
void pyr_fact(float *d_data, float fact1, float *fact2, int32_t n, int32_t nim,
              CarmaDevice *device);

template <class Tin, class Tout>
void roof_abs2(Tout *d_odata, Tin *d_idata, Tout fact, int32_t ns, int32_t nim,
               CarmaDevice *device);

template <class T>
void roof_intensities(T *d_odata, T *d_idata, int32_t *subindx, int32_t *subindy,
                      int32_t ns, int32_t nvalid, int32_t nim, CarmaDevice *device);

template <class T>
void roof_rollmod(T *d_odata, T *d_idata, T *d_mask, float cx, float cy, int32_t np,
                  int32_t ns, CarmaDevice *device);

template <class T>
void roof_fillbin(T *d_odata, T *d_idata, int32_t nrebin, int32_t np, int32_t ns, int32_t nim,
                  CarmaDevice *device);

void copy_imgin_binimg(float *binimg, int32_t *validsubsx, int32_t *validsubsy, int32_t Nb,
                     float *img, int32_t *validx, int32_t *validy, int32_t Nim, int32_t Npix,
                     CarmaDevice *device);

#endif  // _SUTRA_WFS_H_
