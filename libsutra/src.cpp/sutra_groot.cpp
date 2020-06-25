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

//! \file      sutra_groot.cpp
//! \ingroup   libsutra
//! \class     sutra_groot
//! \brief     this class provides the groot features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.4.1
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#include <sutra_groot.h>

sutra_groot::sutra_groot(carma_context *context, int device, int nactus,
                         int nlayers, float gsangle, float *vdt, float *Htheta,
                         float *L0, float *winddir, float *scale, float *pzt2tt,
                         float *TTPfilter, float *Nact, float *xpos,
                         float *ypos, float fc) {
  init_common(context, device, xpos, ypos, nactus, fc);
  this->nactus = nactus;
  this->nlayers = nlayers;
  this->gsangle = gsangle;

  long dims_data1[2] = {1, nlayers};
  this->h_vdt = new carma_host_obj<float>(dims_data1, vdt, MA_PAGELOCK);
  this->h_Htheta = new carma_host_obj<float>(dims_data1, Htheta, MA_PAGELOCK);
  this->h_L0 = new carma_host_obj<float>(dims_data1, L0, MA_PAGELOCK);
  this->h_winddir = new carma_host_obj<float>(dims_data1, winddir, MA_PAGELOCK);
  this->h_scale = new carma_host_obj<float>(dims_data1, scale, MA_PAGELOCK);

  long dims_data2[3] = {2, nactus, nactus};
  this->d_Nact = new carma_obj<float>(context, dims_data2, Nact);
  this->d_Cerr = new carma_obj<float>(context, dims_data2);
  this->d_TTPfilter = new carma_obj<float>(context, dims_data2, TTPfilter);
  dims_data2[1] = 2;
  dims_data2[2] = 2;
  this->d_TT = new carma_obj<float>(context, dims_data2);

  dims_data2[1] = 2;
  dims_data2[2] = nactus;
  this->d_pzt2tt = new carma_obj<float>(context, dims_data2, pzt2tt);

  printf("I am Groot\n");
}

sutra_groot::sutra_groot(carma_context *context, int device, int nssp,
                         float *weights, float scale, float *xpos, float *ypos,
                         float fc, float d, int npts) {
  init_common(context, device, xpos, ypos, nssp, fc);
  this->npts = npts;
  this->d = d;
  this->scale = scale;
  this->nssp = nssp;
  long dims_data1[2] = {1, npts};
  this->h_weights = new carma_host_obj<float>(dims_data1, weights, MA_PAGELOCK);

  long dims_data2[3] = {2, nssp, nssp};
  this->d_CaXX = new carma_obj<float>(context, dims_data2);
  this->d_CaYY = new carma_obj<float>(context, dims_data2);
  printf("I am Groot\n");
}

void sutra_groot::init_common(carma_context *context, int device, float *xpos,
                              float *ypos, int N, float fc) {
  this->current_context = context;
  this->device = device;
  this->current_context->set_activeDevice(device, 1);
  this->fc = fc;

  long dims_data1[2] = {1, N};
  this->d_xpos = new carma_obj<float>(context, dims_data1, xpos);
  this->d_ypos = new carma_obj<float>(context, dims_data1, ypos);

  int Npts = 10000;
  dims_data1[1] = Npts;

  this->d_tab_int_x = new carma_obj<float>(context, dims_data1);
  this->d_tab_int_y = new carma_obj<float>(context, dims_data1);

  tab_u831J0(this->d_tab_int_x->getData(), this->d_tab_int_y->getData(), Npts,
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

sutra_groot::~sutra_groot() {
  this->current_context->set_activeDevice(this->device, 1);
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

int sutra_groot::compute_Cerr() {
  this->current_context->set_activeDevice(this->device, 1);

  carmaSafeCall(cudaMemset(this->d_Cerr->getData(), 0,
                           sizeof(float) * this->d_Cerr->getNbElem()));
  printf("Computing Cerr...\n");

  for (int l = 0; l < this->nlayers; l++) {
    compute_Cerr_layer<float>(
        this->d_Cerr->getData(), this->nactus, this->d_tab_int_x->getData(),
        this->d_tab_int_y->getData(), this->d_xpos->getData(),
        this->d_ypos->getData(), (*this->h_vdt)[l], (*this->h_Htheta)[l],
        (*this->h_L0)[l], this->fc, (*this->h_winddir)[l], this->gsangle,
        (*this->h_scale)[l], this->d_tab_int_y->getNbElem(),
        this->current_context->get_device(device));
  }
  add_transpose<float>(this->d_Cerr->getData(), this->nactus,
                       this->current_context->get_device(device));
  printf("Done\n");

  printf("Applying coupling matrix...\n");
  // Coupling matrix filter
  // carma_magma_potri(this->d_Nact);
  carma_obj<float> d_tmp(this->d_Cerr);
  carma_gemm(this->cublas_handle(), 'n', 'n', this->nactus, this->nactus,
             this->nactus, (float)1., this->d_Nact->getData(), this->nactus,
             this->d_Cerr->getData(), this->nactus, (float)0.0, d_tmp.getData(),
             this->nactus);
  carma_gemm(this->cublas_handle(), 'n', 'n', this->nactus, this->nactus,
             this->nactus, (float)1.0, d_tmp.getData(), this->nactus,
             this->d_Nact->getData(), this->nactus, (float)0.0,
             this->d_Cerr->getData(), this->nactus);
  printf("Done\n");

  // Tip-tilt component
  printf("Computing TT component...\n");
  carma_obj<float> d_tmp2(this->d_pzt2tt);
  carma_gemm(this->cublas_handle(), 'n', 'n', 2, this->nactus, this->nactus,
             (float)1.0, this->d_pzt2tt->getData(), 2, this->d_Cerr->getData(),
             this->nactus, (float)0.0, d_tmp2.getData(), 2);
  carma_gemm(this->cublas_handle(), 'n', 't', 2, 2, this->nactus, (float)1.0,
             d_tmp2.getData(), 2, this->d_pzt2tt->getData(), 2, (float)0.0,
             this->d_TT->getData(), 2);
  printf("Done\n");

  // Filtering TT + piston from Cerr
  printf("Filtering TT + piston from Cerr...\n");
  carma_gemm(this->cublas_handle(), 'n', 'n', this->nactus, this->nactus,
             this->nactus, (float)1.0, this->d_TTPfilter->getData(),
             this->nactus, this->d_Cerr->getData(), this->nactus, (float)0.0,
             d_tmp.getData(), this->nactus);
  carma_gemm(this->cublas_handle(), 'n', 't', this->nactus, this->nactus,
             this->nactus, (float)1.0, d_tmp.getData(), this->nactus,
             this->d_TTPfilter->getData(), this->nactus, (float)0.0,
             this->d_Cerr->getData(), this->nactus);
  printf("Done\n");

  return EXIT_SUCCESS;
}

int sutra_groot::compute_Calias() {
  this->current_context->set_activeDevice(this->device, 1);
  carmaSafeCall(cudaMemset(this->d_CaXX->getData(), 0,
                           sizeof(float) * this->d_CaXX->getNbElem()));
  carmaSafeCall(cudaMemset(this->d_CaYY->getData(), 0,
                           sizeof(float) * this->d_CaYY->getNbElem()));

  float offset;
  for (int k = 0; k < this->npts; k++) {
    offset = k / float(this->npts - 1) * this->d;
    compute_Ca<float>(this->d_CaXX->getData(), this->d_CaYY->getData(),
                      this->nssp, this->d_tab_int_x->getData(),
                      this->d_tab_int_y->getData(), this->d_xpos->getData(),
                      this->d_ypos->getData(), offset, this->d, this->fc,
                      this->scale, this->h_weights->getData()[k],
                      this->d_tab_int_y->getNbElem(),
                      this->current_context->get_device(device));
    if (k > 0) {
      compute_Ca<float>(this->d_CaXX->getData(), this->d_CaYY->getData(),
                        this->nssp, this->d_tab_int_x->getData(),
                        this->d_tab_int_y->getData(), this->d_xpos->getData(),
                        this->d_ypos->getData(), -offset, this->d, this->fc,
                        this->scale, this->h_weights->getData()[k],
                        this->d_tab_int_y->getNbElem(),
                        this->current_context->get_device(device));
    }
  }

  return EXIT_SUCCESS;
}
