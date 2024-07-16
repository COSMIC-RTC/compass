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

//! \file      sutra_groot.cpp
//! \ingroup   libsutra
//! \class     SutraGroot
//! \brief     this class provides the groot features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#include <sutra_groot.hpp>

SutraGroot::SutraGroot(CarmaContext *context, int32_t device, int32_t nactus,
                         int32_t nlayers, float gsangle, float *vdt, float *Htheta,
                         float *L0, float *winddir, float *scale, float *pzt2tt,
                         float *TTPfilter, float *Nact, float *xpos,
                         float *ypos, float fc) {
  init_common(context, device, xpos, ypos, nactus, fc);
  this->nactus = nactus;
  this->nlayers = nlayers;
  this->gsangle = gsangle;

  int64_t dims_data1[2] = {1, nlayers};
  this->h_vdt = new CarmaHostObj<float>(dims_data1, vdt, MA_PAGELOCK);
  this->h_Htheta = new CarmaHostObj<float>(dims_data1, Htheta, MA_PAGELOCK);
  this->h_L0 = new CarmaHostObj<float>(dims_data1, L0, MA_PAGELOCK);
  this->h_winddir = new CarmaHostObj<float>(dims_data1, winddir, MA_PAGELOCK);
  this->h_scale = new CarmaHostObj<float>(dims_data1, scale, MA_PAGELOCK);

  int64_t dims_data2[3] = {2, nactus, nactus};
  this->d_Nact = new CarmaObj<float>(context, dims_data2, Nact);
  this->d_Cerr = new CarmaObj<float>(context, dims_data2);
  this->d_TTPfilter = new CarmaObj<float>(context, dims_data2, TTPfilter);
  dims_data2[1] = 2;
  dims_data2[2] = 2;
  this->d_TT = new CarmaObj<float>(context, dims_data2);

  dims_data2[1] = 2;
  dims_data2[2] = nactus;
  this->d_pzt2tt = new CarmaObj<float>(context, dims_data2, pzt2tt);

  printf("I am Groot\n");
}

SutraGroot::SutraGroot(CarmaContext *context, int32_t device, int32_t nssp,
                         float *weights, float scale, float *xpos, float *ypos,
                         float fc, float d, int32_t npts) {
  init_common(context, device, xpos, ypos, nssp, fc);
  this->npts = npts;
  this->d = d;
  this->scale = scale;
  this->nssp = nssp;
  int64_t dims_data1[2] = {1, npts};
  this->h_weights = new CarmaHostObj<float>(dims_data1, weights, MA_PAGELOCK);

  int64_t dims_data2[3] = {2, nssp, nssp};
  this->d_CaXX = new CarmaObj<float>(context, dims_data2);
  this->d_CaYY = new CarmaObj<float>(context, dims_data2);
  printf("I am Groot\n");
}

void SutraGroot::init_common(CarmaContext *context, int32_t device, float *xpos,
                              float *ypos, int32_t N, float fc) {
  this->current_context = context;
  this->device = device;
  this->current_context->set_active_device(device, 1);
  this->fc = fc;

  int64_t dims_data1[2] = {1, N};
  this->d_xpos = new CarmaObj<float>(context, dims_data1, xpos);
  this->d_ypos = new CarmaObj<float>(context, dims_data1, ypos);

  int32_t Npts = 10000;
  dims_data1[1] = Npts;

  this->d_tab_int_x = new CarmaObj<float>(context, dims_data1);
  this->d_tab_int_y = new CarmaObj<float>(context, dims_data1);

  tab_u831J0(this->d_tab_int_x->get_data(), this->d_tab_int_y->get_data(), Npts,
             this->current_context->get_device(device));
  this->d_CaXX = NULL;
  this->d_CaYY = NULL;
  this->d_Cerr = NULL;
  this->d_TT = NULL;
  this->d_Nact = NULL;
  this->d_pzt2tt = NULL;
  this->d_TTPfilter = NULL;
  this->h_weights = NULL;
  this->h_Htheta = NULL;
  this->h_vdt = NULL;
  this->h_L0 = NULL;
  this->h_winddir = NULL;
  this->h_scale = NULL;
}

SutraGroot::~SutraGroot() {
  this->current_context->set_active_device(this->device, 1);
  if (this->h_vdt != NULL) delete this->h_vdt;
  if (this->h_Htheta != NULL) delete this->h_Htheta;
  if (this->h_L0 != NULL) delete this->h_L0;
  if (this->h_winddir != NULL) delete this->h_winddir;
  if (this->h_scale != NULL) delete this->h_scale;
  if (this->d_xpos != NULL) delete this->d_xpos;
  if (this->d_ypos != NULL) delete this->d_ypos;
  if (this->d_Nact != NULL) delete this->d_Nact;
  if (this->d_TT != NULL) delete this->d_TT;
  if (this->d_Cerr != NULL) delete this->d_Cerr;
  if (this->d_pzt2tt != NULL) delete this->d_pzt2tt;
  if (this->d_TTPfilter != NULL) delete this->d_TTPfilter;
  if (this->d_CaXX != NULL) delete this->d_CaXX;
  if (this->d_CaYY != NULL) delete this->d_CaYY;
  if (this->h_weights != NULL) delete this->h_weights;
  if (this->d_tab_int_x != NULL) delete this->d_tab_int_x;
  if (this->d_tab_int_y != NULL) delete this->d_tab_int_y;
}

int32_t SutraGroot::compute_Cerr() {
  this->current_context->set_active_device(this->device, 1);

  carma_safe_call(cudaMemset(this->d_Cerr->get_data(), 0,
                           sizeof(float) * this->d_Cerr->get_nb_elements()));
  printf("Computing Cerr...\n");

  for (int32_t l = 0; l < this->nlayers; l++) {
    compute_Cerr_layer<float>(
        this->d_Cerr->get_data(), this->nactus, this->d_tab_int_x->get_data(),
        this->d_tab_int_y->get_data(), this->d_xpos->get_data(),
        this->d_ypos->get_data(), (*this->h_vdt)[l], (*this->h_Htheta)[l],
        (*this->h_L0)[l], this->fc, (*this->h_winddir)[l], this->gsangle,
        (*this->h_scale)[l], this->d_tab_int_y->get_nb_elements(),
        this->current_context->get_device(device));
  }
  add_transpose<float>(this->d_Cerr->get_data(), this->nactus,
                       this->current_context->get_device(device));
  printf("Done\n");

  printf("Applying coupling matrix...\n");
  // Coupling matrix filter
  // carma_potr_inv(this->d_Nact);
  CarmaObj<float> d_tmp(this->d_Cerr);
  carma_gemm(this->cublas_handle(), 'n', 'n', this->nactus, this->nactus,
             this->nactus, (float)1., this->d_Nact->get_data(), this->nactus,
             this->d_Cerr->get_data(), this->nactus, (float)0.0, d_tmp.get_data(),
             this->nactus);
  carma_gemm(this->cublas_handle(), 'n', 'n', this->nactus, this->nactus,
             this->nactus, (float)1.0, d_tmp.get_data(), this->nactus,
             this->d_Nact->get_data(), this->nactus, (float)0.0,
             this->d_Cerr->get_data(), this->nactus);
  printf("Done\n");

  // Tip-tilt component
  printf("Computing TT component...\n");
  CarmaObj<float> d_tmp2(this->d_pzt2tt);
  carma_gemm(this->cublas_handle(), 'n', 'n', 2, this->nactus, this->nactus,
             (float)1.0, this->d_pzt2tt->get_data(), 2, this->d_Cerr->get_data(),
             this->nactus, (float)0.0, d_tmp2.get_data(), 2);
  carma_gemm(this->cublas_handle(), 'n', 't', 2, 2, this->nactus, (float)1.0,
             d_tmp2.get_data(), 2, this->d_pzt2tt->get_data(), 2, (float)0.0,
             this->d_TT->get_data(), 2);
  printf("Done\n");

  // Filtering TT + piston from Cerr
  printf("Filtering TT + piston from Cerr...\n");
  carma_gemm(this->cublas_handle(), 'n', 'n', this->nactus, this->nactus,
             this->nactus, (float)1.0, this->d_TTPfilter->get_data(),
             this->nactus, this->d_Cerr->get_data(), this->nactus, (float)0.0,
             d_tmp.get_data(), this->nactus);
  carma_gemm(this->cublas_handle(), 'n', 't', this->nactus, this->nactus,
             this->nactus, (float)1.0, d_tmp.get_data(), this->nactus,
             this->d_TTPfilter->get_data(), this->nactus, (float)0.0,
             this->d_Cerr->get_data(), this->nactus);
  printf("Done\n");

  return EXIT_SUCCESS;
}

int32_t SutraGroot::compute_Calias() {
  this->current_context->set_active_device(this->device, 1);
  carma_safe_call(cudaMemset(this->d_CaXX->get_data(), 0,
                           sizeof(float) * this->d_CaXX->get_nb_elements()));
  carma_safe_call(cudaMemset(this->d_CaYY->get_data(), 0,
                           sizeof(float) * this->d_CaYY->get_nb_elements()));

  float offset;
  for (int32_t k = 0; k < this->npts; k++) {
    offset = k / float(this->npts - 1) * this->d;
    compute_Ca<float>(this->d_CaXX->get_data(), this->d_CaYY->get_data(),
                      this->nssp, this->d_tab_int_x->get_data(),
                      this->d_tab_int_y->get_data(), this->d_xpos->get_data(),
                      this->d_ypos->get_data(), offset, this->d, this->fc,
                      this->scale, this->h_weights->get_data()[k],
                      this->d_tab_int_y->get_nb_elements(),
                      this->current_context->get_device(device));
    if (k > 0) {
      compute_Ca<float>(this->d_CaXX->get_data(), this->d_CaYY->get_data(),
                        this->nssp, this->d_tab_int_x->get_data(),
                        this->d_tab_int_y->get_data(), this->d_xpos->get_data(),
                        this->d_ypos->get_data(), -offset, this->d, this->fc,
                        this->scale, this->h_weights->get_data()[k],
                        this->d_tab_int_y->get_nb_elements(),
                        this->current_context->get_device(device));
    }
  }

  return EXIT_SUCCESS;
}
