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

//! \file      sutra_gamora.hpp
//! \ingroup   libsutra
//! \class     SutraGamora
//! \brief     this class provides the gamora features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#ifndef _SUTRA_GAMORA_H_
#define _SUTRA_GAMORA_H_

#include <carma.hpp>
#include <carma_host_obj.hpp>
#include <carma_obj.hpp>
#include <carma_sparse_obj.hpp>
#include <map>
#include <vector>

class SutraGamora {
 public:
  CarmaContext *current_context;
  int32_t device;

  int32_t nactus;  // number of actuators
  int32_t niter;   // number of iterations
  int32_t nmodes;  // number of modes
  // PSF reconstruction from roket data
  CarmaObj<float> *d_err;
  CarmaObj<cuFloatComplex> *d_amplipup;
  CarmaObj<float> *d_psf;
  CarmaObj<float> *d_phase;
  CarmaObj<int32_t> *d_wherephase;
  CarmaSparseObj<float> *d_IF;
  CarmaObj<float> *d_TT;
  float scale;
  int32_t size;
  int32_t Npts;
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
  SutraGamora(CarmaContext *context, int32_t device, char *type, int32_t nactus,
               int32_t nmodes, int32_t niter, float *IFvalue, int32_t *IFrowind,
               int32_t *IFcolind, int32_t IFnz, float *d_TT, float *pupil, int32_t size,
               int32_t Npts, float scale, float *Btt, float *covmodes);
  ~SutraGamora();

  int32_t psf_rec_roket(float *err);
  int32_t psf_rec_Vii();

  void compute_Dphi_on_mode_k(int32_t k);

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
  std::vector<CarmaObj<int32_t> *> d_wherephase_ngpu;
  std::vector<CarmaObj<cuFloatComplex> *> d_pupfft_ngpu;
  std::vector<CarmaObj<cuFloatComplex> *> d_Dphi_ngpu;
};

int32_t fill_amplipup(cuFloatComplex *amplipup, float *phase, int32_t *wherephase,
                  float scale, int32_t Npts, int32_t nx, int32_t Nx, int32_t puponly,
                  CarmaDevice *device);
int32_t cumulpsf(float *d_odata, cuFloatComplex *d_idata, int32_t N,
             CarmaDevice *device);
int32_t abs2complex(cuFloatComplex *d_odata, cuFloatComplex *d_idata, int32_t N,
                CarmaDevice *device);
int32_t real(float *d_odata, cuFloatComplex *d_idata, int32_t N, CarmaDevice *device);
int32_t fill_mask(float *d_odata, float *d_idata, int32_t N, int32_t norm,
              CarmaDevice *device);
int32_t modulus2(float *d_odata, cuFloatComplex *d_idata, int32_t N,
             CarmaDevice *device);
int32_t pow2(cuFloatComplex *d_odata, cuFloatComplex *d_idata, int32_t N,
         CarmaDevice *device);
int32_t fill_term1(float *d_odata, cuFloatComplex *d_idata,
               cuFloatComplex *d_pupfft, int32_t N, CarmaDevice *device);
int32_t add2Dphi(cuFloatComplex *d_odata, float *d_term1, float *d_term2, float e,
             int32_t N, CarmaDevice *device);
int32_t computeOTFvii(float *d_otfVii, cuFloatComplex *d_Dphi, float *d_otftel,
                  float *d_mask, float scale, int32_t N, CarmaDevice *device);
int32_t ifftscale(cuFloatComplex *d_odata, float scale, int32_t N,
              CarmaDevice *device);

#endif
