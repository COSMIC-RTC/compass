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

//! \file      sutra_gamora.h
//! \ingroup   libsutra
//! \class     sutra_gamora
//! \brief     this class provides the gamora features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.3.2
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#ifndef _SUTRA_GAMORA_H_
#define _SUTRA_GAMORA_H_

#include <carma.h>
#include <carma_host_obj.h>
#include <carma_obj.h>
#include <carma_sparse_obj.h>
#include <map>
#include <vector>

class sutra_gamora {
 public:
  carma_context *current_context;
  int device;

  int nactus;  // number of actuators
  int niter;   // number of iterations
  int nmodes;  // number of modes
  // PSF reconstruction from roket data
  carma_obj<float> *d_err;
  carma_obj<cuFloatComplex> *d_amplipup;
  carma_obj<float> *d_psf;
  carma_obj<float> *d_phase;
  carma_obj<int> *d_wherephase;
  carma_sparse_obj<float> *d_IF;
  carma_obj<float> *d_TT;
  float scale;
  int size;
  int Npts;
  // PSF reconstruction with Vii functions
  carma_obj<float> *d_term1;
  carma_obj<float> *d_term2;
  carma_obj<float> *d_otftel;
  carma_obj<float> *d_otfVii;
  carma_obj<float> *d_mask;
  carma_obj<float> *d_eigenvals;
  carma_obj<float> *d_Btt;
  carma_obj<float> *d_covmodes;
  carma_host_obj<float> *h_eigenvals;
  carma_obj<cuFloatComplex> *d_newmodek;
  carma_obj<cuFloatComplex> *d_Dphi;
  carma_obj<cuFloatComplex> *d_pupfft;

 public:
  sutra_gamora(carma_context *context, int device, char *type, int nactus,
               int nmodes, int niter, float *IFvalue, int *IFrowind,
               int *IFcolind, int IFnz, float *d_TT, float *pupil, int size,
               int Npts, float scale, float *Btt, float *covmodes);
  ~sutra_gamora();

  int psf_rec_roket(float *err);
  int psf_rec_Vii();

  void compute_Dphi_on_mode_k(int k);

 private:
  std::vector<carma_obj<cuFloatComplex> *> d_amplipup_ngpu;
  std::vector<carma_obj<cuFloatComplex> *> d_newmodek_ngpu;
  std::vector<carma_obj<float> *> d_Btt_ngpu;
  std::vector<carma_obj<float> *> d_covmodes_ngpu;
  std::vector<carma_obj<float> *> d_term1_ngpu;
  std::vector<carma_obj<float> *> d_term2_ngpu;
  std::vector<carma_sparse_obj<float> *> d_IF_ngpu;
  std::vector<carma_obj<float> *> d_TT_ngpu;
  std::vector<carma_obj<float> *> d_phase_ngpu;
  std::vector<carma_obj<int> *> d_wherephase_ngpu;
  std::vector<carma_obj<cuFloatComplex> *> d_pupfft_ngpu;
  std::vector<carma_obj<cuFloatComplex> *> d_Dphi_ngpu;
};

int fill_amplipup(cuFloatComplex *amplipup, float *phase, int *wherephase,
                  float scale, int Npts, int nx, int Nx, int puponly,
                  carma_device *device);
int cumulpsf(float *d_odata, cuFloatComplex *d_idata, int N,
             carma_device *device);
int abs2complex(cuFloatComplex *d_odata, cuFloatComplex *d_idata, int N,
                carma_device *device);
int real(float *d_odata, cuFloatComplex *d_idata, int N, carma_device *device);
int fill_mask(float *d_odata, float *d_idata, int N, int norm,
              carma_device *device);
int modulus2(float *d_odata, cuFloatComplex *d_idata, int N,
             carma_device *device);
int pow2(cuFloatComplex *d_odata, cuFloatComplex *d_idata, int N,
         carma_device *device);
int fill_term1(float *d_odata, cuFloatComplex *d_idata,
               cuFloatComplex *d_pupfft, int N, carma_device *device);
int add2Dphi(cuFloatComplex *d_odata, float *d_term1, float *d_term2, float e,
             int N, carma_device *device);
int computeOTFvii(float *d_otfVii, cuFloatComplex *d_Dphi, float *d_otftel,
                  float *d_mask, float scale, int N, carma_device *device);
int ifftscale(cuFloatComplex *d_odata, float scale, int N,
              carma_device *device);

#endif
