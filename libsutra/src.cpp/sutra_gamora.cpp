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

//! \file      sutra_gamora.cpp
//! \ingroup   libsutra
//! \class     SutraGamora
//! \brief     this class provides the gamora features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#include <sutra_gamora.hpp>

using namespace indicators;

SutraGamora::SutraGamora(CarmaContext *context, int32_t device, char *type,
                           int32_t nactus, int32_t nmodes, int32_t niter, float *IFvalue,
                           int32_t *IFrowind, int32_t *IFcolind, int32_t IFnz, float *TT,
                           float *pupil, int32_t size, int32_t Npts, float scale,
                           float *Btt, float *covmodes) {
  this->current_context = context;

  const int32_t ngpu = context->get_ndevice();
  DEBUG_TRACE("using GAMORA with %d GPUs", ngpu);
  if (ngpu == 1)
    this->device = device;
  else {
    int32_t devices[ngpu];
    for (int32_t i = 0; i < ngpu; i++) {
      devices[i] = i;
    }
    this->device = devices[0];
  }

  this->current_context->set_active_device(this->device, 1);

  this->scale = scale;  // 2*pi/lambda, lambda in microns
  this->Npts = Npts;    // number of pupil points
  this->nactus = nactus;
  this->nmodes = nmodes;
  this->niter = niter;
  this->size = size;  // pupdiam in pixels
  this->d_term1 = NULL;
  this->d_term2 = NULL;
  this->d_newmodek = NULL;
  this->d_otftel = NULL;
  this->d_otfVii = NULL;
  this->d_Btt = NULL;
  this->d_covmodes = NULL;
  this->d_mask = NULL;
  this->d_err = NULL;
  this->d_amplipup = NULL;
  this->d_wherephase = NULL;
  this->d_phase = NULL;
  this->d_psf = NULL;
  this->d_IF = NULL;
  this->d_TT = NULL;

  int32_t mradix = 2;
  int32_t fft_size = pow(mradix, (int64_t)(logf(2 * size) / logf(mradix)) + 1);

  int64_t dims_data2[3] = {2, niter, nactus};
  int64_t dims_data1[2] = {1, Npts};

  if (strcmp(type, "roket") == 0) {
    // Command error
    this->d_err = new CarmaObj<float>(this->current_context, dims_data2);
  }

  int32_t *wherephase;
  wherephase = (int32_t *)malloc(Npts * sizeof(int32_t));
  int32_t cpt = 0;
  // Phase point index in spupil
  for (int32_t cc = 0; cc < size * size; cc++) {
    if (pupil[cc] > 0) {
      wherephase[cpt] = cc;
      cpt += 1;
    }
  }
  this->d_wherephase =
      new CarmaObj<int32_t>(this->current_context, dims_data1, wherephase);
  this->d_phase = new CarmaObj<float>(this->current_context, dims_data1);

  // dims_data2[1] = size;
  // dims_data2[2] = size;

  // spupil
  // this->d_pupil = new
  // CarmaObj<float>(this->current_context,dims_data2,pupil);
  // FFT good size
  dims_data2[1] = fft_size;
  dims_data2[2] = fft_size;

  this->d_psf = new CarmaObj<float>(this->current_context, dims_data2);
  this->d_amplipup =
      new CarmaObj<cuFloatComplex>(this->current_context, dims_data2);

  cufftHandle *plan = this->d_amplipup->get_plan();  ///< FFT plan
  carmafft_safe_call(cufftPlan2d(plan, this->d_amplipup->get_dims(1),
                               this->d_amplipup->get_dims(2), CUFFT_C2C));

  dims_data1[1] = IFnz;

  CarmaObj<float> d_val(this->current_context, dims_data1, IFvalue);
  CarmaObj<int32_t> d_row(this->current_context, dims_data1, IFrowind);
  dims_data1[1] = this->nactus - 2 + 1;
  CarmaObj<int32_t> d_col(this->current_context, dims_data1, IFcolind);
  dims_data2[1] = nactus - 2;
  dims_data2[2] = Npts;
  this->d_IF = new CarmaSparseObj<float>(this->current_context, dims_data2,
                                           d_val.get_data(), d_row.get_data(),
                                           d_col.get_data(), IFnz, false);
  dims_data2[1] = 2;
  dims_data2[2] = Npts;
  this->d_TT = new CarmaObj<float>(this->current_context, dims_data2, TT);

  if (strcmp(type, "Vii") == 0) {
    dims_data2[1] = fft_size;
    dims_data2[2] = fft_size;
    this->d_term1 = new CarmaObj<float>(this->current_context, dims_data2);
    this->d_term2 = new CarmaObj<float>(this->current_context, dims_data2);
    this->d_Dphi =
        new CarmaObj<cuFloatComplex>(this->current_context, dims_data2);
    this->d_newmodek =
        new CarmaObj<cuFloatComplex>(this->current_context, dims_data2);
    this->d_pupfft =
        new CarmaObj<cuFloatComplex>(this->current_context, dims_data2);
    this->d_otftel = new CarmaObj<float>(this->current_context, dims_data2);
    this->d_otfVii = new CarmaObj<float>(this->current_context, dims_data2);
    this->d_mask = new CarmaObj<float>(this->current_context, dims_data2);
    dims_data2[1] = nactus;
    dims_data2[2] = nmodes;
    this->d_Btt = new CarmaObj<float>(this->current_context, dims_data2, Btt);
    dims_data2[1] = nmodes;
    this->d_covmodes =
        new CarmaObj<float>(this->current_context, dims_data2, covmodes);
    dims_data1[1] = nmodes;
    this->h_eigenvals = new CarmaHostObj<float>(dims_data1, MA_PAGELOCK);
    this->d_eigenvals = new CarmaObj<float>(this->current_context, dims_data1);
    // FFT plans
    carmafft_safe_call(cufftPlan2d(this->d_pupfft->get_plan(),
                                 this->d_pupfft->get_dims(1),
                                 this->d_pupfft->get_dims(2), CUFFT_C2C));
    carmafft_safe_call(cufftPlan2d(this->d_newmodek->get_plan(),
                                 this->d_newmodek->get_dims(1),
                                 this->d_newmodek->get_dims(2), CUFFT_C2C));
    carmafft_safe_call(cufftPlan2d(this->d_Dphi->get_plan(),
                                 this->d_Dphi->get_dims(1),
                                 this->d_Dphi->get_dims(2), CUFFT_C2C));

    this->d_amplipup_ngpu.push_back(this->d_amplipup);
    this->d_newmodek_ngpu.push_back(this->d_newmodek);
    this->d_Btt_ngpu.push_back(this->d_Btt);
    this->d_covmodes_ngpu.push_back(this->d_covmodes);
    this->d_term1_ngpu.push_back(this->d_term1);
    this->d_term2_ngpu.push_back(this->d_term2);
    this->d_IF_ngpu.push_back(this->d_IF);
    this->d_TT_ngpu.push_back(this->d_TT);
    this->d_phase_ngpu.push_back(this->d_phase);
    this->d_wherephase_ngpu.push_back(this->d_wherephase);
    this->d_pupfft_ngpu.push_back(this->d_pupfft);
    this->d_Dphi_ngpu.push_back(this->d_Dphi);

    for (int32_t d = 1; d < ngpu; d++) {
      current_context->set_active_device(d, 1);
      dims_data2[1] = fft_size;
      dims_data2[2] = fft_size;
      d_amplipup_ngpu.push_back(
          new CarmaObj<cuFloatComplex>(this->current_context, dims_data2));
      cufftHandle *plan = this->d_amplipup_ngpu[d]->get_plan();  ///< FFT plan
      carmafft_safe_call(
          cufftPlan2d(plan, dims_data2[1], dims_data2[2], CUFFT_C2C));
      d_term1_ngpu.push_back(
          new CarmaObj<float>(this->current_context, dims_data2));
      d_term2_ngpu.push_back(
          new CarmaObj<float>(this->current_context, dims_data2));
      d_Dphi_ngpu.push_back(
          new CarmaObj<cuFloatComplex>(this->current_context, dims_data2));
      d_newmodek_ngpu.push_back(
          new CarmaObj<cuFloatComplex>(this->current_context, dims_data2));
      d_pupfft_ngpu.push_back(
          new CarmaObj<cuFloatComplex>(this->current_context, dims_data2));
      carmafft_safe_call(cufftPlan2d(this->d_pupfft_ngpu[d]->get_plan(),
                                   dims_data2[1], dims_data2[2], CUFFT_C2C));
      carmafft_safe_call(cufftPlan2d(this->d_newmodek_ngpu[d]->get_plan(),
                                   dims_data2[1], dims_data2[2], CUFFT_C2C));
      carmafft_safe_call(cufftPlan2d(this->d_Dphi_ngpu[d]->get_plan(),
                                   dims_data2[1], dims_data2[2], CUFFT_C2C));
      dims_data2[1] = nactus;
      dims_data2[2] = nmodes;
      d_Btt_ngpu.push_back(
          new CarmaObj<float>(this->current_context, dims_data2, Btt));
      dims_data2[1] = nmodes;
      d_covmodes_ngpu.push_back(
          new CarmaObj<float>(this->current_context, dims_data2, covmodes));
      dims_data2[1] = nactus - 2;
      dims_data2[2] = Npts;

      dims_data1[1] = IFnz;
      CarmaObj<float> d_val_tmp(this->current_context, dims_data1, IFvalue);
      CarmaObj<int32_t> d_row_tmp(this->current_context, dims_data1, IFrowind);
      dims_data1[1] = this->nactus - 2 + 1;
      CarmaObj<int32_t> d_col_tmp(this->current_context, dims_data1, IFcolind);

      d_IF_ngpu.push_back(new CarmaSparseObj<float>(
          this->current_context, dims_data2, d_val_tmp.get_data(),
          d_row_tmp.get_data(), d_col_tmp.get_data(), IFnz, false));

      dims_data2[1] = 2;
      dims_data2[2] = Npts;
      d_TT_ngpu.push_back(
          new CarmaObj<float>(this->current_context, dims_data2, TT));
      dims_data1[1] = Npts;
      d_wherephase_ngpu.push_back(
          new CarmaObj<int32_t>(this->current_context, dims_data1, wherephase));
      d_phase_ngpu.push_back(
          new CarmaObj<float>(this->current_context, dims_data1));
    }
  }
}

SutraGamora::~SutraGamora() {
  this->current_context->set_active_device(this->device, 1);
  if (this->d_err) delete this->d_err;

  if (this->d_amplipup) {
    for (std::vector<CarmaObj<cuFloatComplex> *>::iterator it =
             this->d_amplipup_ngpu.begin();
         this->d_amplipup_ngpu.end() != it; ++it) {
      current_context->set_active_device((*it)->get_device(), 1);
      delete *it;
    }
  }

  if (this->d_psf) delete this->d_psf;

  if (this->d_pupfft) {
    for (std::vector<CarmaObj<cuFloatComplex> *>::iterator it =
             this->d_pupfft_ngpu.begin();
         this->d_pupfft_ngpu.end() != it; ++it) {
      current_context->set_active_device((*it)->get_device(), 1);
      delete *it;
    }
  }

  if (this->d_phase) {
    for (std::vector<CarmaObj<float> *>::iterator it =
             this->d_phase_ngpu.begin();
         this->d_phase_ngpu.end() != it; ++it) {
      current_context->set_active_device((*it)->get_device(), 1);
      delete *it;
    }
  }

  if (this->d_wherephase) {
    for (std::vector<CarmaObj<int32_t> *>::iterator it =
             this->d_wherephase_ngpu.begin();
         this->d_wherephase_ngpu.end() != it; ++it) {
      current_context->set_active_device((*it)->get_device(), 1);
      delete *it;
    }
  }

  if (this->d_IF) {
    for (std::vector<CarmaSparseObj<float> *>::iterator it =
             this->d_IF_ngpu.begin();
         this->d_IF_ngpu.end() != it; ++it) {
      current_context->set_active_device((*it)->get_device(), 1);
      delete *it;
    }
  }

  if (this->d_term1) {
    for (std::vector<CarmaObj<float> *>::iterator it =
             this->d_term1_ngpu.begin();
         this->d_term1_ngpu.end() != it; ++it) {
      current_context->set_active_device((*it)->get_device(), 1);
      delete *it;
    }

    for (std::vector<CarmaObj<float> *>::iterator it =
             this->d_term2_ngpu.begin();
         this->d_term2_ngpu.end() != it; ++it) {
      current_context->set_active_device((*it)->get_device(), 1);
      delete *it;
    }
    for (std::vector<CarmaObj<cuFloatComplex> *>::iterator it =
             this->d_newmodek_ngpu.begin();
         this->d_newmodek_ngpu.end() != it; ++it) {
      current_context->set_active_device((*it)->get_device(), 1);
      delete *it;
    }
    delete this->d_otftel;
    delete this->d_otfVii;
    delete this->d_mask;
    for (std::vector<CarmaObj<float> *>::iterator it =
             this->d_Btt_ngpu.begin();
         this->d_Btt_ngpu.end() != it; ++it) {
      current_context->set_active_device((*it)->get_device(), 1);
      delete *it;
    }
    for (std::vector<CarmaObj<float> *>::iterator it =
             this->d_covmodes_ngpu.begin();
         this->d_covmodes_ngpu.end() != it; ++it) {
      current_context->set_active_device((*it)->get_device(), 1);
      delete *it;
    }
    delete this->h_eigenvals;
    for (std::vector<CarmaObj<cuFloatComplex> *>::iterator it =
             this->d_Dphi_ngpu.begin();
         this->d_Dphi_ngpu.end() != it; ++it) {
      current_context->set_active_device((*it)->get_device(), 1);
      delete *it;
    }
    for (std::vector<CarmaObj<float> *>::iterator it = this->d_TT_ngpu.begin();
         this->d_TT_ngpu.end() != it; ++it) {
      current_context->set_active_device((*it)->get_device(), 1);
      delete *it;
    }
  }
}

int32_t SutraGamora::psf_rec_roket(float *err) {
  this->current_context->set_active_device(this->device, 1);
  // Get the error command buffer
  this->d_err->host2device(err);
  // set psf to 0
  carma_safe_call(cudaMemset(this->d_psf->get_data(), 0,
                           sizeof(float) * this->d_psf->get_nb_elements()));

  for (int32_t cc = 0; cc < this->niter; cc++) {
    // set amplipup to 0
    carma_safe_call(
        cudaMemset(this->d_amplipup->get_data(), 0,
                   sizeof(cuFloatComplex) * this->d_amplipup->get_nb_elements()));
    // Apply iter #cc on the DM to get residual phase
    carma_gemv(this->current_context->get_cusparse_handle(), 't', 1.0f,
               this->d_IF, this->d_err->get_data_at(cc * this->nactus), 0.0f,
               this->d_phase->get_data());
    carma_gemv(this->current_context->get_cublas_handle(), 't',
               this->d_TT->get_dims(1), this->d_TT->get_dims(2), 1.0f,
               this->d_TT->get_data(), this->d_TT->get_dims(1),
               this->d_err->get_data_at(cc * this->nactus + (this->nactus - 2)),
               1, 1.0f, this->d_phase->get_data(), 1);

    // complex amplitude in the pupil put in the fft support
    fill_amplipup(this->d_amplipup->get_data(), this->d_phase->get_data(),
                  this->d_wherephase->get_data(), this->scale * (-1), this->Npts,
                  this->size, this->d_amplipup->get_dims()[1], 0,
                  this->current_context->get_device(this->device));
    // complex amplitude in the focal plane
    CarmaFFT(this->d_amplipup->get_data(), this->d_amplipup->get_data(), 1,
              *this->d_amplipup->get_plan());
    // take square modulus and add it a the psf LE
    cumulpsf(this->d_psf->get_data(), this->d_amplipup->get_data(),
             this->d_psf->get_dims(1) * this->d_psf->get_dims(2),
             current_context->get_device(device));

    printf("\rComputing and stacking %d PSFs : %d%%", this->niter,
           (cc * 100 / this->niter));
  }
  // rescale because of fft and number of iterations
  this->d_psf->scale(
      1.0f / (this->d_wherephase->get_dims(1) * this->d_wherephase->get_dims(1)),
      1);
  this->d_psf->scale(1.0f / this->niter, 1);
  // DEBUG_TRACE("%d %d",this->Npts*this->Npts*this->niter,this->niter);
  return EXIT_SUCCESS;
}

int32_t SutraGamora::psf_rec_Vii() {
  // Telescope OTF computation and mask
  // Get the pupil
  this->current_context->set_active_device(this->device, 1);

  printf("Computing Telescope OTF and corresponding mask...\n");
  carma_safe_call(
      cudaMemset(this->d_pupfft->get_data(), 0,
                 sizeof(cuFloatComplex) * this->d_pupfft->get_nb_elements()));
  carma_safe_call(cudaMemset(this->d_otftel->get_data(), 0,
                           sizeof(float) * this->d_otftel->get_nb_elements()));

  fill_amplipup(this->d_pupfft->get_data(), this->d_term1->get_data(),
                this->d_wherephase->get_data(), 1.0f, this->Npts, this->size,
                this->d_pupfft->get_dims(1), 1,
                this->current_context->get_device(this->device));
  // compute fft(pupil)

  CarmaFFT(this->d_pupfft->get_data(), this->d_pupfft->get_data(), 1,
            *this->d_pupfft->get_plan());
  // compute pupfft * conjpupfft as abs(pupfft)**2 and store it in
  // the real part of d_amplipup
  abs2complex(this->d_amplipup->get_data(), this->d_pupfft->get_data(),
              this->d_pupfft->get_nb_elements(),
              this->current_context->get_device(this->device));
  // compute ifft(pupfft*conjpupfft)

  CarmaFFT(this->d_amplipup->get_data(), this->d_amplipup->get_data(), -1,
            *this->d_amplipup->get_plan());
  ifftscale(this->d_amplipup->get_data(), 1.0f / this->d_amplipup->get_nb_elements(),
            this->d_amplipup->get_nb_elements(),
            this->current_context->get_device(this->device));
  // Get telescope OTF as real part of amplipup
  real(this->d_otftel->get_data(), this->d_amplipup->get_data(),
       this->d_amplipup->get_nb_elements(),
       this->current_context->get_device(this->device));
  // Deduce the mask from telescope OTF
  fill_mask(this->d_mask->get_data(), this->d_otftel->get_data(),
            this->d_mask->get_nb_elements(), this->Npts,
            this->current_context->get_device(this->device));
  printf("Done\n");

  // OTF|| computation with Vii algorithm
  carma_safe_call(cudaMemset(this->d_Dphi->get_data(), 0,
                           sizeof(cuFloatComplex) * this->d_Dphi->get_nb_elements()));
  carma_safe_call(cudaMemset(this->d_otfVii->get_data(), 0,
                           sizeof(float) * this->d_otfVii->get_nb_elements()));

  // Diagonalisation of covmodes
  carma_syevd<float>(SOLVER_EIG_MODE_VECTOR, this->d_covmodes, this->d_eigenvals);
  d_eigenvals->device2host(h_eigenvals->get_data());

  for (std::vector<CarmaObj<cuFloatComplex> *>::iterator it =
           this->d_pupfft_ngpu.begin();
       this->d_pupfft_ngpu.end() != it; ++it) {
    CarmaObj<cuFloatComplex> *tmp_pupfft = this->d_pupfft;
    if (*it != tmp_pupfft) {
      current_context->set_active_device((*it)->get_device(), 1);
      (*it)->copy_from(tmp_pupfft->get_data(), tmp_pupfft->get_nb_elements());
    }
  }
  for (std::vector<CarmaObj<cuFloatComplex> *>::iterator it =
           this->d_Dphi_ngpu.begin();
       this->d_Dphi_ngpu.end() != it; ++it) {
    if (*it != this->d_Dphi) {
      current_context->set_active_device((*it)->get_device(), 1);
      carma_safe_call(
          cudaMemset((*it)->get_data(), 0, sizeof(float) * (*it)->get_nb_elements()));
    }
  }
  for (std::vector<CarmaObj<float> *>::iterator it =
           this->d_covmodes_ngpu.begin();
       this->d_covmodes_ngpu.end() != it; ++it) {
    CarmaObj<float> *tmp_cov = this->d_covmodes;
    if (*it != tmp_cov) {
      current_context->set_active_device((*it)->get_device(), 1);
      (*it)->copy_from(tmp_cov->get_data(), tmp_cov->get_nb_elements());
    }
  }

  // Loop on the modes to compute OTF from Vii on the fly
  ProgressBar bar{option::BarWidth{50}, option::ForegroundColor{Color::white},
                  option::ShowElapsedTime{true}, option::ShowRemainingTime{true},
                  option::PrefixText{"Computing Vii: "}, option::MaxProgress{this->nmodes}};
  for (int32_t k = 0; k < this->nmodes; k++) {
    compute_Dphi_on_mode_k(k);
    bar.tick();
    // printf("\rComputing OTF with %d Vii :
    // %d%%",this->nmodes,(k*100/this->nmodes));
  }
  this->current_context->set_active_device(this->device, 1);

  CarmaObj<cuFloatComplex> *tmp_vector =
      new CarmaObj<cuFloatComplex>(current_context, this->d_Dphi->get_dims());

  cuFloatComplex alpha;
  alpha.x = 1.0f;
  alpha.y = 0.0f;
  for (std::vector<CarmaObj<cuFloatComplex> *>::iterator it =
           this->d_Dphi_ngpu.begin();
       this->d_Dphi_ngpu.end() != it; ++it) {
    if (*it != d_Dphi) {
      if (current_context->can_p2p(d_Dphi->get_device(), (*it)->get_device())) {
        d_Dphi->axpy(alpha, (*it), 1, 1);
      } else {
        tmp_vector->copy_from((*it)->get_data(), (*it)->get_nb_elements());
        d_Dphi->axpy(alpha, tmp_vector, 1, 1);
      }
    }
  }
  delete tmp_vector;

  // ifft(2*Dphi)
  CarmaFFT(this->d_Dphi->get_data(), this->d_Dphi->get_data(), -1,
            *this->d_Dphi->get_plan());
  ifftscale(this->d_Dphi->get_data(), 1.0f / this->d_Dphi->get_nb_elements(),
            this->d_Dphi->get_nb_elements(),
            this->current_context->get_device(this->device));

  // OTF = exp(-0.5*real(ifft(2*Dphi)) * mask * scale**2 / otftel) * mask
  computeOTFvii(this->d_otfVii->get_data(), this->d_Dphi->get_data(),
                this->d_otftel->get_data(), this->d_mask->get_data(), this->scale,
                this->d_otfVii->get_nb_elements(),
                this->current_context->get_device(this->device));

  return EXIT_SUCCESS;
}

void SutraGamora::compute_Dphi_on_mode_k(int32_t k) {
  this->current_context->set_active_device(this->device, 1);
  int32_t ngpu = d_pupfft_ngpu.size();
  if (ngpu < 2) {
    carma_safe_call(
        cudaMemset(this->d_amplipup->get_data(), 0,
                   sizeof(cuFloatComplex) * this->d_amplipup->get_nb_elements()));
    carma_safe_call(
        cudaMemset(this->d_newmodek->get_data(), 0,
                   sizeof(cuFloatComplex) * this->d_newmodek->get_nb_elements()));
    // Get the new mode shape from d_U, IF and Btt
    carma_gemv<float>(this->current_context->get_cublas_handle(), 'n',
                      this->nactus, this->nmodes, 1.0f, this->d_Btt->get_data(),
                      this->nactus,
                      this->d_covmodes->get_data_at(k * this->nmodes), 1, 0.0f,
                      this->d_term1->get_data(), 1);
    carma_gemv(this->current_context->get_cusparse_handle(), 't', 1.0f,
               this->d_IF, this->d_term1->get_data(), 0.0f,
               this->d_phase->get_data());

    carma_gemv(this->current_context->get_cublas_handle(), 't',
               this->d_TT->get_dims(1), this->d_TT->get_dims(2), 1.0f,
               this->d_TT->get_data(), this->d_TT->get_dims(1),
               this->d_term1->get_data_at(this->nactus - 2), 1, 1.0f,
               this->d_phase->get_data(), 1);
    fill_amplipup(this->d_newmodek->get_data(), this->d_phase->get_data(),
                  this->d_wherephase->get_data(), 1.0f, this->Npts, this->size,
                  this->d_newmodek->get_dims(1), 2,
                  this->current_context->get_device(this->device));
    // Compute term2 = abs(fft(newmode))**2

    CarmaFFT(this->d_newmodek->get_data(), this->d_amplipup->get_data(), 1,
              *this->d_amplipup->get_plan());
    modulus2(this->d_term2->get_data(), this->d_amplipup->get_data(),
             this->d_term2->get_nb_elements(),
             this->current_context->get_device(this->device));
    // Compute term1 = real(fft(newmodek**2)*conjpupfft)
    pow2(this->d_newmodek->get_data(), this->d_newmodek->get_data(),
         this->d_newmodek->get_nb_elements(),
         this->current_context->get_device(this->device));
    CarmaFFT(this->d_newmodek->get_data(), this->d_newmodek->get_data(), 1,
              *this->d_newmodek->get_plan());
    fill_term1(this->d_term1->get_data(), this->d_newmodek->get_data(),
               this->d_pupfft->get_data(), this->d_term1->get_nb_elements(),
               this->current_context->get_device(this->device));
    // Dphi += (term1 - term2) * eigenvals[k]
    add2Dphi(this->d_Dphi->get_data(), this->d_term1->get_data(),
             this->d_term2->get_data(), this->h_eigenvals->get_data()[k],
             this->d_Dphi->get_nb_elements(),
             this->current_context->get_device(this->device));
  } else {
    int32_t cur_device = k % ngpu;
    this->current_context->set_active_device(cur_device, 1);
    carma_safe_call(
        cudaMemset(this->d_amplipup_ngpu[cur_device]->get_data(), 0,
                   sizeof(cuFloatComplex) *
                       this->d_amplipup_ngpu[cur_device]->get_nb_elements()));
    carma_safe_call(
        cudaMemset(this->d_newmodek_ngpu[cur_device]->get_data(), 0,
                   sizeof(cuFloatComplex) *
                       this->d_newmodek_ngpu[cur_device]->get_nb_elements()));
    // Get the new mode shape from d_U, IF and Btt
    carma_gemv<float>(
        this->current_context->get_cublas_handle(), 'n', this->nactus,
        this->nmodes, 1.0f, this->d_Btt_ngpu[cur_device]->get_data(),
        this->nactus,
        this->d_covmodes_ngpu[cur_device]->get_data_at(k * this->nmodes), 1, 0.0f,
        this->d_term1_ngpu[cur_device]->get_data(), 1);

    carma_gemv(this->current_context->get_cusparse_handle(), 't', 1.0f,
               this->d_IF_ngpu[cur_device],
               this->d_term1_ngpu[cur_device]->get_data(), 0.0f,
               this->d_phase_ngpu[cur_device]->get_data());

    carma_gemv(this->current_context->get_cublas_handle(), 't',
               this->d_TT_ngpu[cur_device]->get_dims(1),
               this->d_TT_ngpu[cur_device]->get_dims(2), 1.0f,
               this->d_TT_ngpu[cur_device]->get_data(),
               this->d_TT_ngpu[cur_device]->get_dims(1),
               this->d_term1_ngpu[cur_device]->get_data_at(this->nactus - 2), 1,
               1.0f, this->d_phase_ngpu[cur_device]->get_data(), 1);

    fill_amplipup(this->d_newmodek_ngpu[cur_device]->get_data(),
                  this->d_phase_ngpu[cur_device]->get_data(),
                  this->d_wherephase_ngpu[cur_device]->get_data(), 1.0f,
                  this->Npts, this->size,
                  this->d_newmodek_ngpu[cur_device]->get_dims(1), 2,
                  this->current_context->get_device(cur_device));
    // Compute term2 = abs(fft(newmode))**2

    CarmaFFT(this->d_newmodek_ngpu[cur_device]->get_data(),
              this->d_amplipup_ngpu[cur_device]->get_data(), 1,
              *this->d_amplipup_ngpu[cur_device]->get_plan());

    modulus2(this->d_term2_ngpu[cur_device]->get_data(),
             this->d_amplipup_ngpu[cur_device]->get_data(),
             this->d_term2_ngpu[cur_device]->get_nb_elements(),
             this->current_context->get_device(cur_device));
    // Compute term1 = real(fft(newmodek**2)*conjpupfft)
    pow2(this->d_newmodek_ngpu[cur_device]->get_data(),
         this->d_newmodek_ngpu[cur_device]->get_data(),
         this->d_newmodek_ngpu[cur_device]->get_nb_elements(),
         this->current_context->get_device(cur_device));

    CarmaFFT(this->d_newmodek_ngpu[cur_device]->get_data(),
              this->d_newmodek_ngpu[cur_device]->get_data(), 1,
              *this->d_newmodek_ngpu[cur_device]->get_plan());

    fill_term1(this->d_term1_ngpu[cur_device]->get_data(),
               this->d_newmodek_ngpu[cur_device]->get_data(),
               this->d_pupfft_ngpu[cur_device]->get_data(),
               this->d_term1_ngpu[cur_device]->get_nb_elements(),
               this->current_context->get_device(cur_device));
    // Dphi += (term1 - term2) * eigenvals[k]
    add2Dphi(this->d_Dphi_ngpu[cur_device]->get_data(),
             this->d_term1_ngpu[cur_device]->get_data(),
             this->d_term2_ngpu[cur_device]->get_data(),
             this->h_eigenvals->get_data()[k],
             this->d_Dphi_ngpu[cur_device]->get_nb_elements(),
             this->current_context->get_device(cur_device));
  }
}
