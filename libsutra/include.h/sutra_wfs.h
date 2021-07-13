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

//! \file      sutra_wfs.h
//! \ingroup   libsutra
//! \class     SutraWfs
//! \brief     this class provides the wfs features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.2.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#ifndef _SUTRA_WFS_H_
#define _SUTRA_WFS_H_

#include <carma_utils.h>
#include <sutra_lgs.h>
#include <sutra_phase.h>
#include <sutra_target.h>
#include <sutra_telemetry.h>
#include <sutra_telescope.h>
#include <sutra_utils.h>

#include <map>
#include <vector>
//#include <sutra_slopes.h>

using std::map;
using std::string;

class SutraWfs {
 public:
  int device;
  string type;
  long nxsub;
  long nvalid;
  long npix;
  long nrebin;
  long nfft;
  long ntot;
  long npup;
  long nphase;
  long nmaxhr;
  long nffthr;
  float subapd;
  float nphot;
  float nphot4imat;
  float noise;
  bool lgs;
  bool kernconv;
  bool roket;
  bool is_low_order;

  bool fakecam;
  int max_flux_per_pix;
  int max_pix_value;

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
  CarmaObj<int> *d_hrmap;
  CarmaObj<uint16_t> *d_camimg;
  CarmaObj<float> *d_dark;
  CarmaObj<float> *d_flat;

  CarmaObj<float> *d_slopes;

  CarmaHostObj<float> *image_telemetry;

  SutraSource *d_gs;
  std::vector<CarmaObj<float> *> d_pupil_ngpu;

  CarmaStreams *streams;
  int nstreams;

  CarmaObj<int> *d_phasemap;
  CarmaObj<int> *d_validsubsx;  // nvalid
  CarmaObj<int> *d_validsubsy;  // nvalid

  CarmaContext *current_context;

  /// MPI stuff
  int offset;
  int nvalid_tot;
  int rank;
  int *displ_bincube;
  int *count_bincube;

 public:
  virtual ~SutraWfs(){};

  int wfs_initgs(CarmaObj<float> *d_lgskern,
                 CarmaObj<cuFloatComplex> *d_ftlgskern,
                 map<vector<int>, cufftHandle *> ftlgskern_plans, float xpos,
                 float ypos, float lambda, float mag, float zerop, long size,
                 float noise, long seed, float G, float thetaML, float dx,
                 float dy);
  int set_pupil(float *pupil);
  int set_binimg(float *binimg, int nElem);
  int set_dark(float *dark, int nElem);
  int set_flat(float *flat, int nElem);
  int set_fakecam(bool fakecam);
  int set_max_flux_per_pix(int max_flux_per_pix);
  int set_max_pix_value(int max_pix_value);

  int load_kernels(float *lgskern);
  int sensor_trace(SutraAtmos *yatmos);
  int sensor_trace(SutraDms *ydm, int rst);
  int sensor_trace(SutraAtmos *atmos, SutraDms *ydms);
  int sensor_trace(int rst);
  int slopes_geom(float *slopes, int type = 0);
  int slopes_geom(int type = 0);

  virtual int fill_binimage(int async) = 0;
  virtual int comp_image(bool noise = true) = 0;

  virtual int define_mpi_rank(int rank, int size) = 0;
  virtual int allocate_buffers(
      map<vector<int>, cufftHandle *> campli_plans,
      map<vector<int>, cufftHandle *> fttotim_plans) = 0;
  int set_noise(float noise, long seed);

 protected:
  virtual int comp_generic() = 0;
  SutraWfs(CarmaContext *context, SutraTelescope *d_tel,
            CarmaObj<cuFloatComplex> *d_camplipup,
            CarmaObj<cuFloatComplex> *d_camplifoc,
            CarmaObj<cuFloatComplex> *d_fttotim, string type, long nxsub,
            long nvalid, long npix, long nphase, long nrebin, long nfft,
            long ntot, long npup, float pdiam, float nphotons, float nphot4imat,
            int lgs, bool fakecam, int max_flux_per_pix, int max_pix_value,
            bool is_low_order, bool roket, int device);
};

// General utilities
int fillcamplipup(cuFloatComplex *amplipup, float *phase, float *offset,
                  float *mask, float scale, int *istart, int *jstart,
                  int *ivalid, int *jvalid, int nphase, int npup, int Nfft,
                  int Ntot, CarmaDevice *device, int offset_phase);
int indexfill(cuFloatComplex *d_odata, cuFloatComplex *d_idata, int *indx,
              int nx, int Nx, int N, CarmaDevice *device);
int fillbincube(float *bcube, cuFloatComplex *hrimage, int *indxpix, int Nfft,
                int Npix, int Nrebin, int Nsub, CarmaDevice *device);
int fillbincube_async(CarmaStreams *streams, float *bcube,
                      cuFloatComplex *hrimage, int *indxpix, int Nfft, int Npix,
                      int Nrebin, int Nsub, CarmaDevice *device);
int fillbinimg(float *bimage, float *bcube, int npix, int nsub, int Nsub,
               int *ivalid, int *jvalid, bool add, CarmaDevice *device);
int fillbinimg_async(CarmaStreams *streams, CarmaObj<float> *bimage,
                     CarmaObj<float> *bcube, int npix, int nsub, int Nsub,
                     int *ivalid, int *jvalid, bool add, CarmaDevice *device);
int fillbinimg_async(CarmaHostObj<float> *image_telemetry, float *bimage,
                     float *bcube, int npix, int nsub, int Nsub, int *ivalid,
                     int *jvalid, int nim, bool add, CarmaDevice *device);
int convolve_cube(cuFloatComplex *d_odata, cuFloatComplex *d_idata, int N,
                  int n, CarmaDevice *device);
template <class T>
int digitalize(T *camimg, float *binimg, float *dark, float *flat,
               int max_flux_per_pix, int max_pix_value, int N, CarmaDevice *device);
// CUDA templates
// this is for cog
template <class T>
void subap_reduce(int size, int threads, int blocks, T *d_idata, T *d_odata,
                  CarmaDevice *device);

template <class T>
void subap_reduce_async(int threads, int blocks, CarmaStreams *streams,
                        T *d_idata, T *d_odata);
// this is for tcog
template <class T>
void subap_reduce(int size, int threads, int blocks, T *d_idata, T *d_odata,
                  T thresh, CarmaDevice *device);
// this is for tcog_new
template <class T>
void subap_reduce_new(int size, int threads, int blocks, T *d_idata, T *d_odata,
                      T thresh, CarmaDevice *device);
// this is for wcog
template <class T>
void subap_reduce(int size, int threads, int blocks, T *d_idata, T *d_odata,
                  T *weights, CarmaDevice *device);

template <class T>
void phase_reduce(int threads, int blocks, T *d_idata, T *d_odata, int *indx,
                  T alpha);

template <class T>
void phase_derive(int size, int threads, int blocks, int n, T *d_idata,
                  T *d_odata, int *indx, T *mask, T alpha, float *fluxPerSub);

template <class Tout, class Tin>
void pyr_getpup(Tout *d_odata, Tin *d_idata, Tout *d_offsets, Tin *d_pup,
                int np, float lambda, CarmaDevice *device);

template <class Tout, class Tin>
void pyr_getpup(Tout *d_odata, Tin *d_idata, Tin *d_pup, int np, int N,
                float lambda, float cx, float cy, CarmaDevice *device);

template <class T>
void pyr_rollmod(T *d_odata, T *d_idata, T *d_mask, float cx, float cy, int np,
                 int ns, CarmaDevice *device);

template <class T>
void pyr_fillbin(T *d_odata, T *d_idata, int nrebin, int np, int ns, int nim,
                 CarmaDevice *device);

template <class T>
int pyr_fillbinimg(T *bimage, const T *bcube, const int nxsub, const bool add,
                   CarmaDevice *device);

template <class T>
int pyr_fillbinimg(T *oimage, const T *image, const int n, const int N,
                   const int rebin, const bool add, CarmaDevice *device);

template <class Tin, class Tout>
void pyr_abs2(Tout *d_odata, Tin *d_idata, Tout fact, int ns, int nim,
              CarmaDevice *device);

template <class Tout, class Tin>
void pyr_submask(Tout *d_odata, Tin *d_mask, int n, CarmaDevice *device);

template <class T>
void pyr_submaskpyr(T *d_odata, T *d_mask, int n, CarmaDevice *device);

template <class T>
void pyr_submaskpyr(T *d_odata, float *d_mask, int n, CarmaDevice *device);

template <class Tout, class Tin>
void pyr_abs(Tout *d_odata, Tin *d_idata, int ns, int nim,
             CarmaDevice *device);

template <class Tout, class Tin>
void pyr_submask3d(Tout *d_odata, Tin *d_mask, int n, int nim,
                   CarmaDevice *device);

template <class T>
void pyr_intensities(T *d_odata, T *d_idata, int *subindx, int *subindy, int ns,
                     int nvalid, int nim, CarmaDevice *device);

template <class T>
void pyr_intensities(T *d_odata, T *d_idata, int *subindx, int *subindy, int ns,
                     int nvalid, CarmaDevice *device);

template <class T>
void pyr_fact(T *d_data, T fact, int n, int nim, CarmaDevice *device);

void pyr_fact(cuFloatComplex *d_data, float fact, int n, int nim,
              CarmaDevice *device);
void pyr_fact(float *d_data, float fact1, float *fact2, int n, int nim,
              CarmaDevice *device);

template <class Tin, class Tout>
void roof_abs2(Tout *d_odata, Tin *d_idata, Tout fact, int ns, int nim,
               CarmaDevice *device);

template <class T>
void roof_intensities(T *d_odata, T *d_idata, int *subindx, int *subindy,
                      int ns, int nvalid, int nim, CarmaDevice *device);

template <class T>
void roof_rollmod(T *d_odata, T *d_idata, T *d_mask, float cx, float cy, int np,
                  int ns, CarmaDevice *device);

template <class T>
void roof_fillbin(T *d_odata, T *d_idata, int nrebin, int np, int ns, int nim,
                  CarmaDevice *device);

void copy_imgin_binimg(float *binimg, int *validsubsx, int *validsubsy, int Nb,
                     float *img, int *validx, int *validy, int Nim, int Npix,
                     CarmaDevice *device);

#endif  // _SUTRA_WFS_H_
