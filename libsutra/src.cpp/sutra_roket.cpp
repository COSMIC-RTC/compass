// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
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

//! \file      sutra_roket.cpp
//! \ingroup   libsutra
//! \class     SutraRoket
//! \brief     this class provides the roket features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.2.1
//! \date      2022/01/24
//! \copyright GNU Lesser General Public License

#include <sutra_roket.h>

SutraRoket::SutraRoket(CarmaContext *context, int device, SutraRtc *rtc,
                         SutraSensors *sensors, SutraTarget *target,
                         SutraDms *dms, SutraTelescope *tel, SutraAtmos *atm,
                         int loopcontroller, int geocontroller, int nactus,
                         int nmodes, int nfilt, int niter, float *Btt, float *P,
                         float *gRD, float *RD) {
  // context
  this->current_context = context;
  this->device = device;
  this->current_context->set_active_device(device, 1);
  // sutra objects to supervise
  this->rtc = rtc;
  this->sensors = sensors;
  this->target = target;
  this->dms = dms;
  this->atm = atm;
  this->tel = tel;
  this->loopcontrol = 0L;
  this->geocontrol = 0L;
  // simple initialisation
  this->niter = niter;
  this->nfilt = nfilt;
  this->nactus = nactus;
  this->nmodes = nmodes;
  this->iterk = 0;
  this->loopcontroller = loopcontroller;
  this->geocontroller = geocontroller;
  if (this->rtc->d_control[this->loopcontroller]->get_type().compare("ls") ==
      0) {
    this->loopcontrol = dynamic_cast<sutra_controller_ls *>(
        this->rtc->d_control[this->loopcontroller]);
  }
  if (this->rtc->d_control[this->geocontroller]->get_type().compare("geo") ==
      0) {
    this->geocontrol = dynamic_cast<sutra_controller_geo *>(
        this->rtc->d_control[this->geocontroller]);
  }
  this->gain = loopcontrol->gain;
  this->fitting = 0.;
  this->nslopes = this->loopcontrol->d_centroids->get_dims(1);
  // contributors buffers initialsations
  long dims_data2[3] = {2, this->niter, this->nactus};
  this->d_noise = new CarmaObj<float>(this->current_context, dims_data2);
  this->d_nonlinear = new CarmaObj<float>(this->current_context, dims_data2);
  this->d_tomo = new CarmaObj<float>(this->current_context, dims_data2);
  this->d_filtered = new CarmaObj<float>(this->current_context, dims_data2);
  this->d_alias = new CarmaObj<float>(this->current_context, dims_data2);
  this->d_bandwidth = new CarmaObj<float>(this->current_context, dims_data2);
  this->d_commanded = new CarmaObj<float>(this->current_context, dims_data2);

  // Matrix initialisation
  dims_data2[1] = this->nmodes;
  this->d_P = new CarmaObj<float>(this->current_context, dims_data2, P);
  dims_data2[1] = this->nactus;
  this->d_gRD = new CarmaObj<float>(this->current_context, dims_data2, gRD);
  this->d_RD = new CarmaObj<float>(this->current_context, dims_data2, RD);
  dims_data2[2] = this->nmodes;
  this->d_Btt = new CarmaObj<float>(this->current_context, dims_data2, Btt);
  dims_data2[1] = this->nactus;
  dims_data2[2] = this->nactus;
  this->d_covv = new CarmaObj<float>(this->current_context, dims_data2);
  dims_data2[1] = this->nslopes;
  dims_data2[2] = this->nslopes;
  this->d_covm = new CarmaObj<float>(this->current_context, dims_data2);

  carma_safe_call(cudaMemset(this->d_covv->get_data(), 0,
                           sizeof(float) * this->d_covv->get_nb_elements()));
  carma_safe_call(cudaMemset(this->d_covm->get_data(), 0,
                           sizeof(float) * this->d_covm->get_nb_elements()));
  carma_safe_call(cudaMemset(this->d_noise->get_data(), 0,
                           sizeof(float) * this->d_noise->get_nb_elements()));
  carma_safe_call(cudaMemset(this->d_nonlinear->get_data(), 0,
                           sizeof(float) * this->d_nonlinear->get_nb_elements()));
  carma_safe_call(cudaMemset(this->d_tomo->get_data(), 0,
                           sizeof(float) * this->d_tomo->get_nb_elements()));
  carma_safe_call(cudaMemset(this->d_filtered->get_data(), 0,
                           sizeof(float) * this->d_filtered->get_nb_elements()));
  carma_safe_call(cudaMemset(this->d_alias->get_data(), 0,
                           sizeof(float) * this->d_alias->get_nb_elements()));
  carma_safe_call(cudaMemset(this->d_bandwidth->get_data(), 0,
                           sizeof(float) * this->d_bandwidth->get_nb_elements()));
  carma_safe_call(cudaMemset(this->d_commanded->get_data(), 0,
                           sizeof(float) * this->d_commanded->get_nb_elements()));
  // Residual error buffer initialsations
  long dims_data1[2] = {1, this->nactus};
  this->d_fullErr = new CarmaObj<float>(this->current_context, dims_data1);
  this->d_err1 = new CarmaObj<float>(this->current_context, dims_data1);
  this->d_err2 = new CarmaObj<float>(this->current_context, dims_data1);
  this->d_bkup_com = new CarmaObj<float>(this->current_context, dims_data1);
  this->d_tmpdiff = new CarmaObj<float>(this->current_context, dims_data1);

  // Additional buffers initialsations
  dims_data1[1] = this->nmodes;
  this->d_modes = new CarmaObj<float>(this->current_context, dims_data1);
  this->d_filtmodes = new CarmaObj<float>(this->current_context, dims_data1);

  // Target screen backup
  this->d_bkup_screen = new CarmaObj<float>(
      this->current_context,
      this->target->d_targets[0]->d_phase->d_screen->get_dims());
  // PSF fitting
  this->d_psfortho = new CarmaObj<float>(
      this->current_context, this->target->d_targets[0]->d_image_se->get_dims());
  // this->d_psfse = new
  // CarmaObj<float>(this->current_context,this->target->d_targets[0]->d_image_se->get_dims());
  // carma_safe_call(cudaMemset(this->d_psfse->get_data(), 0, sizeof(float) *
  // this->d_psfse->get_nb_elements()));
}

SutraRoket::~SutraRoket() {
  this->current_context->set_active_device(this->device, 1);
  if (this->d_covv) delete this->d_covv;
  if (this->d_covm) delete this->d_covm;
  if (this->d_P) delete this->d_P;
  if (this->d_Btt) delete this->d_Btt;
  if (this->d_noise) delete this->d_noise;
  if (this->d_nonlinear) delete this->d_nonlinear;
  if (this->d_tomo) delete this->d_tomo;
  if (this->d_filtered) delete this->d_filtered;
  if (this->d_alias) delete this->d_alias;
  if (this->d_bandwidth) delete this->d_bandwidth;
  if (this->d_fullErr) delete this->d_fullErr;
  if (this->d_err1) delete this->d_err1;
  if (this->d_err2) delete this->d_err2;
  if (this->d_bkup_screen) delete this->d_bkup_screen;
  if (this->d_bkup_com) delete this->d_bkup_com;
  if (this->d_commanded) delete this->d_commanded;
  if (this->d_modes) delete this->d_modes;
  if (this->d_filtmodes) delete this->d_filtmodes;
  if (this->d_gRD) delete this->d_gRD;
  if (this->d_RD) delete this->d_RD;
  if (this->d_tmpdiff) delete this->d_tmpdiff;
  if (this->d_psfortho) delete this->d_psfortho;
  //    if(this->d_psfse)
  //        delete this->d_psfse;
}

int SutraRoket::save_loop_state() {
  this->current_context->set_active_device(this->device, 1);
  this->d_fullErr->copy_from(this->loopcontrol->d_err->get_data(), this->nactus);
  this->d_bkup_com->copy_from(this->loopcontrol->d_com->get_data(), this->nactus);
  this->d_bkup_screen->copy_from(
      this->target->d_targets[0]->d_phase->d_screen->get_data(),
      this->target->d_targets[0]->d_phase->d_screen->get_nb_elements());
  // Stack covariance matrices
  carma_ger(this->current_context->get_cublas_handle(), this->nactus,
            this->nactus, 1.0f / this->niter,
            this->loopcontrol->d_com->get_data(), 1,
            this->loopcontrol->d_com->get_data(), 1, this->d_covv->get_data(),
            this->nactus);
  carma_ger(this->current_context->get_cublas_handle(), this->nslopes,
            this->nslopes, 1.0f / this->niter,
            this->loopcontrol->d_centroids->get_data(), 1,
            this->loopcontrol->d_centroids->get_data(), 1,
            this->d_covm->get_data(), this->nslopes);

  return EXIT_SUCCESS;
}

int SutraRoket::restore_loop_state() {
  this->current_context->set_active_device(this->device, 1);
  this->d_bkup_com->copy_into(this->loopcontrol->d_com->get_data(), this->nactus);
  this->d_bkup_screen->copy_into(
      this->target->d_targets[0]->d_phase->d_screen->get_data(),
      this->target->d_targets[0]->d_phase->d_screen->get_nb_elements());

  return EXIT_SUCCESS;
}

int SutraRoket::apply_loop_filter(CarmaObj<float> *d_odata,
                                   CarmaObj<float> *d_idata1,
                                   CarmaObj<float> *d_idata2, float gain,
                                   int k) {
  this->current_context->set_active_device(this->device, 1);
  this->d_tmpdiff->copy_from(d_idata1->get_data(), this->nactus);
  this->d_tmpdiff->axpy(-1.0f, d_idata2, 1, 1);  // err1-err2
  this->d_tmpdiff->copy_into(d_odata->get_data_at((k + 1) * this->nactus),
                            this->nactus);  // odata[k+1,:] = err1 - err2
  carma_gemv<float>(this->current_context->get_cublas_handle(), 'n',
                    this->nactus, this->nactus, 1.0f, this->d_gRD->get_data(),
                    this->nactus, d_odata->get_data_at(k * this->nactus), 1, gain,
                    d_odata->get_data_at((k + 1) * this->nactus),
                    1);  // odata[k+1,:] = gRD*odata[k,:] + g*(err1-err2)

  return EXIT_SUCCESS;
}

int SutraRoket::compute_breakdown() {
  this->current_context->set_active_device(this->device, 1);
  save_loop_state();
  // Noise
  this->rtc->do_centroids(this->loopcontroller, false);
  this->rtc->do_control(this->loopcontroller);
  this->d_err1->copy_from(this->loopcontrol->d_err->get_data(), this->nactus);
  apply_loop_filter(this->d_noise, this->d_fullErr, this->d_err1, this->gain,
                    this->iterk);
  // Non-linearity
  this->rtc->do_centroids_geom(this->loopcontroller);
  this->rtc->do_control(this->loopcontroller);
  this->d_err2->copy_from(this->loopcontrol->d_err->get_data(), this->nactus);
  apply_loop_filter(this->d_nonlinear, this->d_err1, this->d_err2, this->gain,
                    this->iterk);
  // Aliasing on GS direction
  this->geocontrol->comp_dphi(this->sensors->d_wfs[0]->d_gs, true);
  this->rtc->do_control(this->geocontroller);
  this->rtc->apply_control(this->geocontroller);
  this->sensors->d_wfs[0]->sensor_trace(this->dms, 0);
  this->sensors->d_wfs[0]->d_gs->d_phase->d_screen->axpy(
      1.0, this->tel->d_phase_ab_M1_m, 1, 1);
  this->sensors->d_wfs[0]->comp_image();
  this->rtc->do_centroids(this->loopcontroller, false);
  this->rtc->do_control(this->loopcontroller);
  this->d_err1->copy_from(this->loopcontrol->d_err->get_data(), this->nactus);
  this->d_err2->copy_from(this->d_tmpdiff->get_data(), this->nactus);
  apply_loop_filter(this->d_alias, this->d_err1, this->d_err2, this->gain,
                    this->iterk);
  // Wavefront
  this->target->d_targets[0]->raytrace(this->atm);
  this->target->d_targets[0]->d_phase->d_screen->axpy(
      1.0, this->tel->d_phase_ab_M1_m, 1, 1);
  this->current_context->set_active_device(this->geocontrol->device, 1);
  this->geocontrol->comp_dphi(this->target->d_targets[0], false);

  this->rtc->do_control(this->geocontroller);
  this->d_err1->copy_from(this->geocontrol->d_com->get_data(), this->nactus);
  // Fitting
  this->rtc->apply_control(this->geocontroller);

  this->target->d_targets[0]->raytrace(this->dms, 0, 0);
  this->fitting += this->target->d_targets[0]->phase_var / this->niter;
  this->target->d_targets[0]->comp_image(0, false);
  // this->d_psfse->copy_from(this->target->d_targets[0]->d_image_se->get_data(),this->d_psfse->get_nb_elements());
  this->d_psfortho->axpy(1.0f / this->niter,
                         this->target->d_targets[0]->d_image_se, 1, 1);
  // Filtered modes
  carma_gemv<float>(this->current_context->get_cublas_handle(), 'n',
                    this->nmodes, this->nactus, 1.0f, this->d_P->get_data(),
                    this->nmodes, this->d_err1->get_data(), 1, 0.f,
                    this->d_modes->get_data(), 1);
  separate_modes(this->d_modes->get_data(), this->d_filtmodes->get_data(),
                 this->nmodes, this->nfilt,
                 this->current_context->get_device(this->device));
  carma_gemv<float>(this->current_context->get_cublas_handle(), 'n',
                    this->nactus, this->nmodes, 1.0f, this->d_Btt->get_data(),
                    this->nactus, this->d_filtmodes->get_data(), 1, 0.f,
                    this->d_filtered->get_data_at(this->iterk * this->nactus), 1);
  // Commanded modes
  carma_gemv<float>(this->current_context->get_cublas_handle(), 'n',
                    this->nactus, this->nmodes, 1.0f, this->d_Btt->get_data(),
                    this->nactus, this->d_modes->get_data(), 1, 0.f,
                    this->d_commanded->get_data_at(this->iterk * this->nactus),
                    1);
  // Bandwidth
  if (this->iterk > 0) {
    this->d_err1->copy_from(
        this->d_commanded->get_data_at(this->iterk * this->nactus), this->nactus);
    this->d_err2->copy_from(
        this->d_commanded->get_data_at((this->iterk - 1) * this->nactus),
        this->nactus);
    apply_loop_filter(this->d_bandwidth, this->d_err1, this->d_err2, -1.0f,
                      this->iterk - 1);
  } else {
    carma_safe_call(
        cudaMemset(this->d_err2->get_data(), 0, sizeof(float) * this->nactus));
    this->d_err2->axpy(-1.0f, this->d_err1, 1, 1);
    this->d_err2->copy_into(this->d_bandwidth->get_data(), this->nactus);
  }
  // tomography
  this->sensors->d_wfs[0]->sensor_trace(this->atm);
  this->geocontrol->comp_dphi(this->sensors->d_wfs[0]->d_gs, true);
  this->rtc->do_control(this->geocontroller);
  this->d_err1->copy_from(this->geocontrol->d_com->get_data(), this->nactus);
  carma_gemv<float>(this->current_context->get_cublas_handle(), 'n',
                    this->nmodes, this->nactus, 1.0f, this->d_P->get_data(),
                    this->nmodes, this->d_err1->get_data(), 1, 0.f,
                    this->d_modes->get_data(), 1);
  separate_modes(this->d_modes->get_data(), this->d_filtmodes->get_data(),
                 this->nmodes, this->nfilt,
                 this->current_context->get_device(this->device));
  carma_gemv<float>(this->current_context->get_cublas_handle(), 'n',
                    this->nactus, this->nmodes, 1.0f, this->d_Btt->get_data(),
                    this->nactus, this->d_modes->get_data(), 1, 0.f,
                    this->d_err2->get_data(), 1);
  this->d_err1->copy_from(
      this->d_commanded->get_data_at(this->iterk * this->nactus), this->nactus);
  this->d_err1->axpy(-1.0f, this->d_err2, 1, 1);
  carma_gemv<float>(this->current_context->get_cublas_handle(), 'n',
                    this->nactus, this->nactus, 1.0f, this->d_RD->get_data(),
                    this->nactus, this->d_err1->get_data(), 1, 0.f,
                    this->d_err2->get_data(), 1);
  carma_safe_call(
      cudaMemset(this->d_err1->get_data(), 0, sizeof(float) * this->nactus));
  apply_loop_filter(this->d_tomo, this->d_err2, this->d_err1, this->gain,
                    this->iterk);

  restore_loop_state();
  this->iterk++;

  return EXIT_SUCCESS;
}
