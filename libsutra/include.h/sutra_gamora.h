// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_gamora.h
//! \ingroup   libsutra
//! \class     SutraGamora
//! \brief     this class provides the gamora features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#ifndef _SUTRA_GAMORA_H_
#define _SUTRA_GAMORA_H_

#include <carma.h>
#include <carma_host_obj.h>
#include <carma_obj.h>
#include <carma_sparse_obj.h>
#include <map>
#include <vector>

class SutraGamora {
 public:
  CarmaContext *current_context;
  int device;

  int nactus;  // number of actuators
  int niter;   // number of iterations
  int nmodes;  // number of modes
  // PSF reconstruction from roket data
  CarmaObj<float> *d_err;
  CarmaObj<cuFloatComplex> *d_amplipup;
  CarmaObj<float> *d_psf;
  CarmaObj<float> *d_phase;
  CarmaObj<int> *d_wherephase;
  CarmaSparseObj<float> *d_IF;
  CarmaObj<float> *d_TT;
  float scale;
  int size;
  int Npts;
  // PSF reconstruction with Vii functions
  CarmaObj<float> *d_term1;
  CarmaObj<float> *d_term2;
  CarmaObj<float> *d_otftel;
  CarmaObj<float> *d_otfVii;
  CarmaObj<float> *d_mask;
  CarmaObj<float> *d_eigenvals;
  CarmaObj<float> *d_Btt;
  CarmaObj<float> *d_covmodes;
  CarmaHostObj<float> *h_eigenvals;
  CarmaObj<cuFloatComplex> *d_newmodek;
  CarmaObj<cuFloatComplex> *d_Dphi;
  CarmaObj<cuFloatComplex> *d_pupfft;

 public:
  SutraGamora(CarmaContext *context, int device, char *type, int nactus,
               int nmodes, int niter, float *IFvalue, int *IFrowind,
               int *IFcolind, int IFnz, float *d_TT, float *pupil, int size,
               int Npts, float scale, float *Btt, float *covmodes);
  ~SutraGamora();

  int psf_rec_roket(float *err);
  int psf_rec_Vii();

  void compute_Dphi_on_mode_k(int k);

 private:
  std::vector<CarmaObj<cuFloatComplex> *> d_amplipup_ngpu;
  std::vector<CarmaObj<cuFloatComplex> *> d_newmodek_ngpu;
  std::vector<CarmaObj<float> *> d_Btt_ngpu;
  std::vector<CarmaObj<float> *> d_covmodes_ngpu;
  std::vector<CarmaObj<float> *> d_term1_ngpu;
  std::vector<CarmaObj<float> *> d_term2_ngpu;
  std::vector<CarmaSparseObj<float> *> d_IF_ngpu;
  std::vector<CarmaObj<float> *> d_TT_ngpu;
  std::vector<CarmaObj<float> *> d_phase_ngpu;
  std::vector<CarmaObj<int> *> d_wherephase_ngpu;
  std::vector<CarmaObj<cuFloatComplex> *> d_pupfft_ngpu;
  std::vector<CarmaObj<cuFloatComplex> *> d_Dphi_ngpu;
};

int fill_amplipup(cuFloatComplex *amplipup, float *phase, int *wherephase,
                  float scale, int Npts, int nx, int Nx, int puponly,
                  CarmaDevice *device);
int cumulpsf(float *d_odata, cuFloatComplex *d_idata, int N,
             CarmaDevice *device);
int abs2complex(cuFloatComplex *d_odata, cuFloatComplex *d_idata, int N,
                CarmaDevice *device);
int real(float *d_odata, cuFloatComplex *d_idata, int N, CarmaDevice *device);
int fill_mask(float *d_odata, float *d_idata, int N, int norm,
              CarmaDevice *device);
int modulus2(float *d_odata, cuFloatComplex *d_idata, int N,
             CarmaDevice *device);
int pow2(cuFloatComplex *d_odata, cuFloatComplex *d_idata, int N,
         CarmaDevice *device);
int fill_term1(float *d_odata, cuFloatComplex *d_idata,
               cuFloatComplex *d_pupfft, int N, CarmaDevice *device);
int add2Dphi(cuFloatComplex *d_odata, float *d_term1, float *d_term2, float e,
             int N, CarmaDevice *device);
int computeOTFvii(float *d_otfVii, cuFloatComplex *d_Dphi, float *d_otftel,
                  float *d_mask, float scale, int N, CarmaDevice *device);
int ifftscale(cuFloatComplex *d_odata, float scale, int N,
              CarmaDevice *device);

#endif
