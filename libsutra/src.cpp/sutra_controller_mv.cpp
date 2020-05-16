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

//! \file      sutra_controller_mv.cpp
//! \ingroup   libsutra
//! \class     sutra_controller_mv
//! \brief     this class provides the controller_mv features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.0.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#include <carma_magma.h>
#include <sutra_controller_mv.h>
#include <sutra_controller_utils.h>
#include <string>

template <typename Tcomp, typename Tout>
sutra_controller_mv<Tcomp, Tout>::sutra_controller_mv(
    CarmaContext *context, long nvalid_, long nslope_, long nactu_,
    float delay, SutraDms *dms, int *idx_dms, int ndm, int *idx_centro, int ncentro)
    : SutraController<Tcomp, Tout>(context, nvalid_, nslope_, nactu_, delay,
                                    dms, idx_dms, ndm, idx_centro, ncentro) {
  this->gain = 0.0f;

  //  this->nstreams = 1; //nvalid/10;
  //  while (this->nactu() % this->nstreams != 0)
  //    nstreams--;
  //  std::cerr << "controller uses " << nstreams << " streams" << std::endl;
  //  streams = new CarmaStreams(nstreams);
  long dims_data2[3];
  dims_data2[0] = 2;
  dims_data2[1] = this->nslope();
  dims_data2[2] = this->nactu();
  this->d_imat = new CarmaObj<Tcomp>(this->current_context, dims_data2);
  dims_data2[1] = this->nactu();
  dims_data2[2] = this->nslope();
  this->d_cmat = new CarmaObj<Tcomp>(this->current_context, dims_data2);
  // dims_data2[1] = dims_data2[2] = this->nactu();
  // d_U = new CarmaObj<Tcomp>(this->current_context, dims_data2);
  this->d_cenbuff = 0L;
  if ((int)delay > 0) {
    dims_data2[1] = this->nslope();
    dims_data2[2] = (int)delay + 1;
    this->d_cenbuff = new CarmaObj<Tcomp>(this->current_context, dims_data2);
  }
  dims_data2[1] = this->nslope();
  dims_data2[2] = this->nslope();
  this->d_Cmm = new CarmaObj<Tcomp>(this->current_context, dims_data2);

  long dims_data1[2];
  dims_data1[0] = 1;

  // dims_data1[1] = this->nslope() < this->nactu() ? this->nslope() :
  // this->nactu();
  this->h_eigenvals = 0L;
  this->h_Cmmeigenvals = 0L;

  dims_data1[1] = this->nslope();
  this->d_noisemat = new CarmaObj<Tcomp>(this->current_context, dims_data1);
  this->d_olmeas = new CarmaObj<Tcomp>(this->current_context, dims_data1);
  // this->d_compbuff = new CarmaObj<Tcomp>(this->current_context, dims_data1);
  this->d_compbuff2 = new CarmaObj<Tcomp>(this->current_context, dims_data1);
  dims_data1[1] = this->nactu();
  this->d_compbuff = new CarmaObj<Tcomp>(this->current_context, dims_data1);
  this->d_err = new CarmaObj<Tcomp>(this->current_context, dims_data1);
  this->d_gain = new CarmaObj<Tcomp>(this->current_context, dims_data1);

  // Florian features
  dims_data2[1] = this->nactu();  // 64564;
  dims_data2[2] = this->nactu();
  this->d_covmat = new CarmaObj<Tcomp>(this->current_context, dims_data2);
  dims_data2[1] = this->nactu();
  dims_data2[2] = this->nactu();
  this->d_KLbasis = new CarmaObj<Tcomp>(this->current_context, dims_data2);
  this->d_Cphim = 0L;

  cublas_handle = this->current_context->get_cublas_handle();
  // carma_checkCublasStatus(cublasCreate(&(this->cublas_handle)));
}

template <typename Tcomp, typename Tout>
sutra_controller_mv<Tcomp, Tout>::~sutra_controller_mv() {
  this->current_context->set_active_device(this->device, 1);
  // delete this->d_U;

  delete this->d_imat;
  delete this->d_cmat;
  delete this->d_gain;
  delete this->d_Cmm;
  if (this->d_Cphim != 0L) delete this->d_Cphim;
  delete this->d_olmeas;
  delete this->d_compbuff;
  delete this->d_compbuff2;
  delete this->d_noisemat;

  delete this->h_eigenvals;
  delete this->h_Cmmeigenvals;

  if (this->delay > 0) delete this->d_cenbuff;
  delete this->d_err;
  // Florian features
  delete this->d_covmat;
  delete this->d_KLbasis;

  // carma_checkCublasStatus(cublasDestroy(this->cublas_handle));

  // delete this->current_context;
}

template <typename Tcomp, typename Tout>
string sutra_controller_mv<Tcomp, Tout>::get_type() {
  return "mv";
}

template <typename Tcomp, typename Tout>
int sutra_controller_mv<Tcomp, Tout>::set_modal_gains(Tcomp *mgain) {
  this->current_context->set_active_device(this->device, 1);
  this->d_gain->host2device(mgain);
  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int sutra_controller_mv<Tcomp, Tout>::set_cmat(Tcomp *cmat) {
  this->current_context->set_active_device(this->device, 1);
  this->d_cmat->host2device(cmat);
  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int sutra_controller_mv<Tcomp, Tout>::set_imat(Tcomp *imat) {
  this->current_context->set_active_device(this->device, 1);
  this->d_imat->host2device(imat);
  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int sutra_controller_mv<Tcomp, Tout>::load_noisemat(Tcomp *noise) {
  this->current_context->set_active_device(this->device, 1);
  this->d_noisemat->host2device(noise);
  return EXIT_SUCCESS;
}
// Florian features
template <typename Tcomp, typename Tout>
int sutra_controller_mv<Tcomp, Tout>::compute_Cmm(
    SutraAtmos *atmos, SutraSensors *sensors, double *L0, double *cn2,
    double *alphaX, double *alphaY, double diamTel, double cobs) {
  this->current_context->set_active_device(this->device, 1);

  struct gtomo_struct g_tomo;
  init_tomo_gpu_gb(&g_tomo, atmos, sensors, diamTel, cobs);
  update_tomo_sys_gpu_gb(&g_tomo, sensors, alphaX, alphaY);
  update_tomo_atm_gpu_gb(&g_tomo, sensors, atmos, L0, cn2, alphaX, alphaY);
  matcov_gpu_4(this->d_Cmm->get_data(), this->nslope(), this->nslope(), 0, 0,
               this->nslope(), &g_tomo, atmos, sensors, alphaX, alphaY);
  free_tomo_gpu_gb(&g_tomo);

  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int sutra_controller_mv<Tcomp, Tout>::compute_Cphim(
    SutraAtmos *atmos, SutraSensors *sensors, SutraDms *dms, double *L0,
    double *cn2, double *alphaX, double *alphaY, double *X, double *Y,
    double *xactu, double *yactu, double diamTel, double *k2, long *NlayerDm,
    long *indLayerDm, double FoV, double *pitch, double *alt_dm) {
  this->current_context->set_active_device(this->device, 1);

  // Find number of actuators without TTcomp DM
  int Nactu = 0;
  vector<SutraDm *>::iterator p;
  p = dms->d_dms.begin();
  while (p != dms->d_dms.end()) {
    SutraDm *dm = *p;
    if (dm->type != "tt") {
      Nactu += dm->nactus;
    }
    p++;
  }

  long dims_data2[3] = {2, Nactu, this->nslope()};
  if (this->d_Cphim != 0L) delete this->d_Cphim;
  this->d_Cphim = new CarmaObj<Tcomp>(this->current_context, dims_data2);

  struct cphim_struct cphim_struct;

  // Compute Cphim matrix
  init_cphim_struct(&cphim_struct, atmos, sensors, dms, diamTel);
  update_cphim_sys(&cphim_struct, sensors, alphaX, alphaY, xactu, yactu, X, Y,
                   NlayerDm, indLayerDm, alt_dm, pitch, k2, FoV);
  update_cphim_atm(&cphim_struct, sensors, atmos, L0, cn2, alphaX, alphaY);
  CPHIM(this->d_Cphim->get_data(), Nactu, this->nslope(), 0, 0, Nactu,
        &cphim_struct, atmos, sensors, alphaX, alphaY,
        this->current_context->get_device(this->device));
  free_cphim_struct(&cphim_struct);

  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int sutra_controller_mv<Tcomp, Tout>::filter_cphim(Tcomp *F, Tcomp *Nact) {
  // Piston filter
  piston_filt_cphim(this->d_Cphim, F);
  // Init and inverse the coupling matrix
  long dims_data2[3] = {2, this->d_Cphim->get_dims()[1],
                        this->d_Cphim->get_dims()[1]};
  CarmaObj<Tcomp> *d_Nact =
      new CarmaObj<Tcomp>(this->current_context, dims_data2);
  dims_data2[2] = this->nslope();
  CarmaObj<Tcomp> *d_tmp =
      new CarmaObj<Tcomp>(this->current_context, dims_data2);
  d_Nact->host2device(Nact);
  carma_magma_potri(d_Nact);
  carma_gemm(cublas_handle, 'n', 'n', this->d_Cphim->get_dims()[1],
             this->nslope(), this->d_Cphim->get_dims()[1], 1.0f,
             d_Nact->get_data(), this->d_Cphim->get_dims()[1],
             this->d_Cphim->get_data(), this->d_Cphim->get_dims()[1], 0.0f,
             d_tmp->get_data(), this->d_Cphim->get_dims()[1]);
  this->d_Cphim->copy(d_tmp, 1, 1);

  delete d_Nact;
  delete d_tmp;

  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int sutra_controller_mv<Tcomp, Tout>::do_covmat(SutraDm *ydm, char *method,
                                                int *indx_pup, long dim,
                                                Tcomp *xpos, Tcomp *ypos,
                                                long Nkl, Tcomp norm,
                                                Tcomp ampli) {
  this->current_context->set_active_device(this->device, 1);
  long dims_data[3];
  dims_data[0] = 2;
  dims_data[1] = Nkl;
  dims_data[2] = Nkl;
  CarmaObj<Tcomp> *d_statcov =
      new CarmaObj<Tcomp>(this->current_context, dims_data);
  long dims_data2[2];
  dims_data2[0] = 1;
  dims_data2[1] = this->nactu();
  CarmaObj<Tcomp> *d_KLcov =
      new CarmaObj<Tcomp>(this->current_context, dims_data2);

  dims_data2[1] = dim;
  CarmaObj<int> *d_indx =
      new CarmaObj<int>(this->current_context, dims_data2);

  // Compute the statistic matrix from actuators positions & Kolmogorov
  // statistic
  CarmaObj<Tcomp> *d_xpos =
      new CarmaObj<Tcomp>(this->current_context, dims_data);
  CarmaObj<Tcomp> *d_ypos =
      new CarmaObj<Tcomp>(this->current_context, dims_data);
  d_xpos->host2device(xpos);
  d_ypos->host2device(ypos);
  do_statmat(d_statcov->get_data(), Nkl, d_xpos->get_data(), d_ypos->get_data(),
             norm, this->current_context->get_device(this->device));
  // Compute and apply piston filter
  this->piston_filt(d_statcov);

  if (ydm->type == "pzt") {
    this->d_covmat->copy(d_statcov, 1, 1);
    delete d_statcov;
    delete d_KLcov;
    delete d_indx;
    delete d_xpos;
    delete d_ypos;

    return EXIT_SUCCESS;
    /*
     dims_data[1] = dim;
     dims_data[2] = this->nactu();
     CarmaObj<Tcomp> *d_IF = new CarmaObj<Tcomp>(this->current_context,
     dims_data); dims_data[1] = this->nactu(); dims_data[2] = this->nactu();
     CarmaObj<Tcomp> *d_geocov = new CarmaObj<Tcomp>(current_context,
     dims_data);

     // Get influence functions of the DM (to be CUsparsed)
     d_indx->host2device(indx_pup);
     ydm->get_IF(d_IF->get_data(),d_indx->get_data(),dim);

     // Compute geometric matrix (to be CUsparsed)
     this->do_geomat(d_geocov,d_IF,dim,ampli);

     delete d_IF;

     // Double diagonalisation to obtain KL basis on actuators
     this->DDiago(d_statcov,d_geocov);

     delete d_geocov;
     */
  }
  // Computation of covariance matrix
  // 1. Computation of covariance matrix in KL basis
  CarmaHostObj<Tcomp> *h_eigenvals =
      new CarmaHostObj<Tcomp>(dims_data2, MA_PAGELOCK);

  carma_magma_syevd<Tcomp>('N', d_statcov, h_eigenvals);

  if (ydm->type == "kl") {
    dims_data2[1] = this->nactu();
    CarmaHostObj<Tcomp> h_KLcov(dims_data2, MA_PAGELOCK);
    if (strcmp(method, "inv") == 0) {
      for (int i = 0; i < this->nactu(); i++) {
        h_KLcov[i] = -1. / (h_eigenvals->get_data())[i];
      }
      if (Nkl == this->nactu()) {
        h_KLcov[this->nactu() - 1] = 0.;
      }
      d_KLcov->host2device(h_KLcov.get_data());
      add_md(this->d_covmat->get_data(), this->d_covmat->get_data(),
             d_KLcov->get_data(), this->nactu(),
             this->current_context->get_device(this->device));
    }
    if (strcmp(method, "n") == 0) {
      for (int i = 0; i < this->nactu(); i++) {
        h_KLcov[i] = -(h_eigenvals->get_data())[i];
      }
      d_KLcov->host2device(h_KLcov.get_data());
      add_md(this->d_covmat->get_data(), this->d_covmat->get_data(),
             d_KLcov->get_data(), this->nactu(),
             this->current_context->get_device(this->device));
    }
    delete h_KLcov;
  }
  if (ydm->type == "pzt") {
    if (strcmp(method, "inv") == 0) {
      // Inversion of the KL covariance matrix
      CarmaHostObj<Tcomp> *h_eigenvals_inv =
          new CarmaHostObj<Tcomp>(dims_data2, MA_PAGELOCK);
      for (int i = 0; i < this->nactu(); i++) {
        // Valeurs propres négatives.... A voir et inverser ordre si valeurs
        // propres positives
        h_eigenvals_inv->get_data()[i] = -1. / h_eigenvals->get_data()[i];
      }

      h_eigenvals_inv->get_data()[this->nactu() - 1] = 0.;
      d_KLcov->host2device(h_eigenvals_inv->get_data());

      // 2. Inversion of the KL basis
      dims_data2[0] = 1;
      dims_data2[1] = this->nactu();
      CarmaObj<Tcomp> *d_eigen =
          new CarmaObj<Tcomp>(this->current_context, dims_data2);
      dims_data[1] = this->nactu();
      dims_data[2] = this->nactu();
      CarmaObj<Tcomp> *d_tmp =
          new CarmaObj<Tcomp>(this->current_context, dims_data);
      CarmaObj<Tcomp> *d_Ukl =
          new CarmaObj<Tcomp>(this->current_context, dims_data);
      CarmaObj<Tcomp> *d_Vkl =
          new CarmaObj<Tcomp>(this->current_context, dims_data);
      CarmaHostObj<Tcomp> *h_KL =
          new CarmaHostObj<Tcomp>(dims_data, MA_PAGELOCK);
      CarmaHostObj<Tcomp> *h_U =
          new CarmaHostObj<Tcomp>(dims_data, MA_PAGELOCK);
      CarmaHostObj<Tcomp> *h_Vt =
          new CarmaHostObj<Tcomp>(dims_data, MA_PAGELOCK);

      h_KL->cpy_obj(this->d_KLbasis, cudaMemcpyDeviceToHost);

      carma_magma_svd_cpu<Tcomp>(h_KL, h_eigenvals, h_U, h_Vt);

      d_Ukl->host2device(h_Vt->get_data());
      d_Vkl->host2device(h_U->get_data());

      for (int i = 0; i < this->nactu(); i++) {
        h_eigenvals_inv->get_data()[i] = 1. / h_eigenvals->get_data()[i];
      }

      d_eigen->host2device(h_eigenvals_inv->get_data());

      carma_dgmm(cublas_handle, CUBLAS_SIDE_LEFT, this->nactu(), this->nactu(),
                 d_Vkl->get_data(), this->nactu(), d_eigen->get_data(), 1,
                 d_tmp->get_data(), this->nactu());
      carma_gemm(cublas_handle, 't', 't', this->nactu(), this->nactu(),
                 this->nactu(), 1.0f, d_tmp->get_data(), this->nactu(),
                 d_Ukl->get_data(), this->nactu(), 0.0f,
                 this->d_KLbasis->get_data(), this->nactu());

      // 3. Projection of KL covariance matrix in the DM basis

      carma_dgmm(cublas_handle, CUBLAS_SIDE_LEFT, this->nactu(), this->nactu(),
                 d_KLbasis->get_data(), this->nactu(), d_KLcov->get_data(), 1,
                 d_tmp->get_data(), this->nactu());
      carma_gemm(cublas_handle, 't', 'n', this->nactu(), this->nactu(),
                 this->nactu(), 1.0f, d_tmp->get_data(), this->nactu(),
                 this->d_KLbasis->get_data(), this->nactu(), 0.0f,
                 this->d_covmat->get_data(), this->nactu());

      delete d_eigen;
      delete d_tmp;
      delete d_Ukl;
      delete d_Vkl;
      delete h_KL;
      delete h_U;
      delete h_Vt;
      delete h_eigenvals_inv;
    } else if (strcmp(method, "n") == 0) {
      dims_data[1] = this->nactu();
      dims_data[2] = this->nactu();
      CarmaObj<Tcomp> *d_tmp =
          new CarmaObj<Tcomp>(this->current_context, dims_data);
      for (int i = 0; i < this->nactu(); i++) {
        // Valeurs propres négatives.... A voir et inverser ordre si valeurs
        // propres positives
        h_eigenvals->get_data()[i] = -h_eigenvals->get_data()[i];
        std::cout << h_eigenvals->get_data()[i] << std::endl;
      }
      h_eigenvals->get_data()[this->nactu() - 1] = 0.;
      d_KLcov->host2device(h_eigenvals->get_data());

      carma_dgmm(cublas_handle, CUBLAS_SIDE_RIGHT, this->nactu(), this->nactu(),
                 d_KLbasis->get_data(), this->nactu(), d_KLcov->get_data(), 1,
                 d_tmp->get_data(), this->nactu());
      carma_gemm(cublas_handle, 'n', 't', this->nactu(), this->nactu(),
                 this->nactu(), 1.0f, d_tmp->get_data(), this->nactu(),
                 this->d_KLbasis->get_data(), this->nactu(), 0.0f,
                 this->d_covmat->get_data(), this->nactu());

      delete d_tmp;
    }
  }

  delete d_statcov;
  delete h_eigenvals;
  delete d_KLcov;
  delete d_indx;
  delete d_xpos;
  delete d_ypos;

  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int sutra_controller_mv<Tcomp, Tout>::do_geomat(CarmaObj<Tcomp> *d_geocov,
                                                CarmaObj<Tcomp> *d_IF,
                                                long n_pts, Tcomp ampli) {
  this->current_context->set_active_device(this->device, 1);
  carma_gemm(cublas_handle, 't', 'n', this->nactu(), this->nactu(), n_pts, 1.0f,
             d_IF->get_data(), n_pts, d_IF->get_data(), n_pts, 0.0f,
             d_geocov->get_data(), this->nactu());
  d_geocov->scale(ampli, 1);
  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int sutra_controller_mv<Tcomp, Tout>::piston_filt(CarmaObj<Tcomp> *d_statcov) {
  this->current_context->set_active_device(this->device, 1);
  long Nmod = d_statcov->get_dims()[1];
  long dims_data[3];
  dims_data[0] = 2;
  dims_data[1] = Nmod;
  dims_data[2] = Nmod;
  CarmaObj<Tcomp> *d_F =
      new CarmaObj<Tcomp>(this->current_context, dims_data);
  CarmaObj<Tcomp> *d_tmp =
      new CarmaObj<Tcomp>(this->current_context, dims_data);

  int N = d_statcov->get_dims()[1] * d_statcov->get_dims()[1];
  fill_filtmat(d_F->get_data(), Nmod, N,
               this->current_context->get_device(this->device));

  carma_gemm(cublas_handle, 'n', 'n', Nmod, Nmod, Nmod, 1.0f, d_F->get_data(),
             Nmod, d_statcov->get_data(), Nmod, 0.0f, d_tmp->get_data(), Nmod);
  carma_gemm(cublas_handle, 'n', 'n', Nmod, Nmod, Nmod, 1.0f, d_tmp->get_data(),
             Nmod, d_F->get_data(), Nmod, 0.0f, d_statcov->get_data(), Nmod);

  delete d_tmp;
  delete d_F;

  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int sutra_controller_mv<Tcomp, Tout>::piston_filt_cphim(
    CarmaObj<Tcomp> *d_cphim, Tcomp *F) {
  this->current_context->set_active_device(this->device, 1);

  long Nmod = d_cphim->get_dims()[1];
  long dims_data[3];
  dims_data[0] = 2;
  dims_data[1] = Nmod;
  dims_data[2] = Nmod;
  CarmaObj<Tcomp> *d_F =
      new CarmaObj<Tcomp>(this->current_context, dims_data);
  dims_data[2] = this->nslope();
  CarmaObj<Tcomp> *d_tmp =
      new CarmaObj<Tcomp>(this->current_context, dims_data);

  d_F->host2device(F);

  carma_gemm(cublas_handle, 'n', 'n', Nmod, this->nslope(), Nmod, 1.0f,
             d_F->get_data(), Nmod, d_cphim->get_data(), Nmod, 0.0f,
             d_tmp->get_data(), Nmod);
  d_cphim->copy(d_tmp, 1, 1);

  delete d_tmp;
  delete d_F;

  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int sutra_controller_mv<Tcomp, Tout>::invgen(CarmaObj<Tcomp> *d_mat,
                                             Tcomp cond, int job) {
  this->current_context->set_active_device(this->device, 1);
  const long dims_data[3] = {2, d_mat->get_dims()[1], d_mat->get_dims()[2]};
  CarmaObj<Tcomp> *d_U =
      new CarmaObj<Tcomp>(this->current_context, dims_data);
  CarmaObj<Tcomp> *d_tmp =
      new CarmaObj<Tcomp>(this->current_context, dims_data);
  int i;

  const long dims_data2[2] = {1, d_mat->get_dims()[1]};
  CarmaObj<Tcomp> *d_eigenvals_inv =
      new CarmaObj<Tcomp>(this->current_context, dims_data2);
  CarmaHostObj<Tcomp> *h_eigenvals =
      new CarmaHostObj<Tcomp>(dims_data2, MA_PAGELOCK);
  CarmaHostObj<Tcomp> *h_eigenvals_inv =
      new CarmaHostObj<Tcomp>(dims_data2, MA_PAGELOCK);

  d_U->copy(d_mat, 1, 1);
  carma_magma_syevd<Tcomp>('V', d_U, h_eigenvals);
  // syevd_f('V',d_U,h_eigenvals);
  if (job == 1) {  // Conditionnement
    Tcomp maxe = h_eigenvals->get_data()[d_mat->get_dims()[1] - 1];

    for (i = 0; i < d_mat->get_dims()[1]; i++) {
      if (h_eigenvals->get_data()[i] < maxe / cond)
        h_eigenvals_inv->get_data()[i] = 0.;
      else
        h_eigenvals_inv->get_data()[i] = 1. / h_eigenvals->get_data()[i];
    }
  }
  if (job == 0) {  // Filtre #cond modes
    for (i = 0; i < d_mat->get_dims()[1]; i++) {
      if (i < cond)
        h_eigenvals_inv->get_data()[i] = 0.;
      else
        h_eigenvals_inv->get_data()[i] = 1. / h_eigenvals->get_data()[i];
    }
  }

  d_eigenvals_inv->host2device(h_eigenvals_inv->get_data());

  carma_dgmm(cublas_handle, CUBLAS_SIDE_RIGHT, d_mat->get_dims()[1],
             d_mat->get_dims()[2], d_U->get_data(), d_mat->get_dims()[1],
             d_eigenvals_inv->get_data(), 1, d_tmp->get_data(),
             d_mat->get_dims()[1]);
  carma_gemm<Tcomp>(cublas_handle, 'n', 't', d_mat->get_dims()[1],
                    d_mat->get_dims()[1], d_mat->get_dims()[2], 1.0f,
                    d_tmp->get_data(), d_mat->get_dims()[1], d_U->get_data(),
                    d_mat->get_dims()[1], 0.0f, d_mat->get_data(),
                    d_mat->get_dims()[1]);

  delete d_U;
  delete d_tmp;
  delete d_eigenvals_inv;
  delete h_eigenvals;
  delete h_eigenvals_inv;

  return EXIT_SUCCESS;
}
template <typename Tcomp, typename Tout>
int sutra_controller_mv<Tcomp, Tout>::invgen(CarmaObj<Tcomp> *d_mat,
                                             CarmaHostObj<Tcomp> *h_eigen,
                                             Tcomp cond) {
  this->current_context->set_active_device(this->device, 1);
  const long dims_data[3] = {2, d_mat->get_dims()[1], d_mat->get_dims()[2]};
  CarmaObj<Tcomp> *d_U =
      new CarmaObj<Tcomp>(this->current_context, dims_data);
  CarmaObj<Tcomp> *d_tmp =
      new CarmaObj<Tcomp>(this->current_context, dims_data);

  const long dims_data2[2] = {1, d_mat->get_dims()[1]};
  CarmaObj<Tcomp> *d_eigenvals_inv =
      new CarmaObj<Tcomp>(this->current_context, dims_data2);
  CarmaHostObj<Tcomp> *h_eigenvals_inv =
      new CarmaHostObj<Tcomp>(dims_data2, MA_PAGELOCK);

  d_U->copy(d_mat, 1, 1);

  carma_magma_syevd<Tcomp>('V', d_U, h_eigen);

  // syevd_f('V',d_U,h_eigen);
  // Conditionnement
  Tcomp maxe = h_eigen->get_data()[d_mat->get_dims()[1] - 3];
  int cpt = 0;
  for (int i = 0; i < d_mat->get_dims()[1]; i++) {
    if (h_eigen->get_data()[i] < maxe / cond) {
      h_eigenvals_inv->get_data()[i] = 0.;
      cpt++;
    } else
      h_eigenvals_inv->get_data()[i] = 1. / h_eigen->get_data()[i];
  }

  d_eigenvals_inv->host2device(h_eigenvals_inv->get_data());

  carma_dgmm(cublas_handle, CUBLAS_SIDE_RIGHT, d_mat->get_dims()[1],
             d_mat->get_dims()[2], d_U->get_data(), d_mat->get_dims()[1],
             d_eigenvals_inv->get_data(), 1, d_tmp->get_data(),
             d_mat->get_dims()[1]);
  carma_gemm<Tcomp>(cublas_handle, 'n', 't', d_mat->get_dims()[1],
                    d_mat->get_dims()[1], d_mat->get_dims()[2], 1.0f,
                    d_tmp->get_data(), d_mat->get_dims()[1], d_U->get_data(),
                    d_mat->get_dims()[1], 0.0f, d_mat->get_data(),
                    d_mat->get_dims()[1]);

  std::cout << "Inversion done with " << cpt << " modes filtered" << std::endl;

  delete d_U;
  delete d_tmp;
  delete d_eigenvals_inv;
  delete h_eigenvals_inv;

  return EXIT_SUCCESS;
}
template <typename Tcomp, typename Tout>
int sutra_controller_mv<Tcomp, Tout>::invgen_cpu(CarmaObj<Tcomp> *d_mat,
                                                 CarmaHostObj<Tcomp> *h_eigen,
                                                 Tcomp cond) {
  this->current_context->set_active_device(this->device, 1);
  const long dims_data[3] = {2, d_mat->get_dims()[1], d_mat->get_dims()[2]};
  CarmaObj<Tcomp> *d_U =
      new CarmaObj<Tcomp>(this->current_context, dims_data);
  CarmaHostObj<Tcomp> *h_U =
      new CarmaHostObj<Tcomp>(dims_data, MA_PAGELOCK);
  CarmaHostObj<Tcomp> *h_V =
      new CarmaHostObj<Tcomp>(dims_data, MA_PAGELOCK);
  CarmaHostObj<Tcomp> *h_mat =
      new CarmaHostObj<Tcomp>(dims_data, MA_PAGELOCK);
  CarmaObj<Tcomp> *d_tmp =
      new CarmaObj<Tcomp>(this->current_context, dims_data);

  const long dims_data2[2] = {1, d_mat->get_dims()[1]};
  CarmaObj<Tcomp> *d_eigenvals_inv =
      new CarmaObj<Tcomp>(this->current_context, dims_data2);
  CarmaHostObj<Tcomp> *h_eigenvals_inv =
      new CarmaHostObj<Tcomp>(dims_data2, MA_PAGELOCK);

  d_mat->device2host(h_mat->get_data());

  carma_magma_svd_cpu<Tcomp>(h_mat, h_eigen, h_U, h_V);
  d_U->host2device(h_V->get_data());
  // syevd_f('V',d_U,h_eigen);
  // Conditionnement
  Tcomp maxe = h_eigen->get_data()[2];
  int cpt = 0;
  for (int i = 0; i < d_mat->get_dims()[1]; i++) {
    if (h_eigen->get_data()[i] < maxe / cond) {
      h_eigenvals_inv->get_data()[i] = 0.;
      cpt++;
    } else
      h_eigenvals_inv->get_data()[i] = 1. / h_eigen->get_data()[i];
  }

  d_eigenvals_inv->host2device(h_eigenvals_inv->get_data());

  carma_dgmm(cublas_handle, CUBLAS_SIDE_RIGHT, d_mat->get_dims()[1],
             d_mat->get_dims()[2], d_U->get_data(), d_mat->get_dims()[1],
             d_eigenvals_inv->get_data(), 1, d_tmp->get_data(),
             d_mat->get_dims()[1]);
  carma_gemm<Tcomp>(cublas_handle, 'n', 't', d_mat->get_dims()[1],
                    d_mat->get_dims()[1], d_mat->get_dims()[2], 1.0f,
                    d_tmp->get_data(), d_mat->get_dims()[1], d_U->get_data(),
                    d_mat->get_dims()[1], 0.0f, d_mat->get_data(),
                    d_mat->get_dims()[1]);

  std::cout << "Inversion done with " << cpt << " modes filtered" << std::endl;

  delete d_U;
  delete d_tmp;
  delete d_eigenvals_inv;
  delete h_eigenvals_inv;
  delete h_U;
  delete h_V;
  delete h_mat;

  return EXIT_SUCCESS;
}

/*
template<typename T>
 int sutra_controller_mv<Tcomp, Tout>::do_statmat(Tcomp *statcov, Tcomp *xpos,
Tcomp *ypos){ int dim_x = sizeof(xpos)/sizeof(xpos[0]); int ind; for (i=0 ;
i<dim_x ; i++){ for(j=0 ; j<dim_x ; j++){ ind = i*dim_x + j; statcov[ind] = 6.88
* pow(sqrt(pow(x[i]-x[j],2) + pow(y[i]-y[j],2)),5./3);
 }
 }

 return EXIT_SUCCESS;
 }
 */
template <typename Tcomp, typename Tout>
int sutra_controller_mv<Tcomp, Tout>::DDiago(CarmaObj<Tcomp> *d_statcov,
                                             CarmaObj<Tcomp> *d_geocov) {
  this->current_context->set_active_device(this->device, 1);
  const long dims_data[3] = {2, this->nactu(), this->nactu()};
  CarmaObj<Tcomp> *d_M1 =
      new CarmaObj<Tcomp>(this->current_context, dims_data);
  CarmaObj<Tcomp> *d_tmp =
      new CarmaObj<Tcomp>(this->current_context, dims_data);
  CarmaObj<Tcomp> *d_tmp2 =
      new CarmaObj<Tcomp>(this->current_context, dims_data);

  const long dims_data2[2] = {1, this->nactu()};
  CarmaObj<Tcomp> *d_eigenvals =
      new CarmaObj<Tcomp>(this->current_context, dims_data2);
  CarmaObj<Tcomp> *d_eigenvals_sqrt =
      new CarmaObj<Tcomp>(this->current_context, dims_data2);
  CarmaObj<Tcomp> *d_eigenvals_inv =
      new CarmaObj<Tcomp>(this->current_context, dims_data2);
  CarmaHostObj<Tcomp> *h_eigenvals =
      new CarmaHostObj<Tcomp>(dims_data2, MA_PAGELOCK);
  CarmaHostObj<Tcomp> *h_eigenvals_inv =
      new CarmaHostObj<Tcomp>(dims_data2, MA_PAGELOCK);
  CarmaHostObj<Tcomp> *h_eigenvals_sqrt =
      new CarmaHostObj<Tcomp>(dims_data2, MA_PAGELOCK);

  // 1. SVdec(geocov,U) --> Ut * geocov * U = D²
  carma_magma_syevd<Tcomp>('V', d_geocov, h_eigenvals);

  d_eigenvals->host2device(h_eigenvals->get_data());
  for (int i = 0; i < this->nactu(); i++) {
    h_eigenvals_sqrt->get_data()[i] =
        sqrt(h_eigenvals->get_data()[i]);  // D = sqrt(D²)
    h_eigenvals_inv->get_data()[i] =
        1. / sqrt(h_eigenvals->get_data()[i]);  // D⁻¹ = 1/sqrt(D²)
  }
  d_eigenvals_sqrt->host2device(h_eigenvals_sqrt->get_data());
  d_eigenvals_inv->host2device(h_eigenvals_inv->get_data());

  // 2. M⁻¹ = sqrt(eigenvals) * Ut : here, we have transpose(M⁻¹)

  carma_dgmm<Tcomp>(cublas_handle, CUBLAS_SIDE_RIGHT, this->nactu(),
                    this->nactu(), d_geocov->get_data(), this->nactu(),
                    d_eigenvals_sqrt->get_data(), 1, d_M1->get_data(),
                    this->nactu());

  // 3. C' = M⁻¹ * statcov * M⁻¹t
  carma_gemm<Tcomp>(cublas_handle, 't', 'n', this->nactu(), this->nactu(),
                    this->nactu(), 1.0f, d_M1->get_data(), this->nactu(),
                    d_statcov->get_data(), this->nactu(), 0.0f, d_tmp->get_data(),
                    this->nactu());

  carma_gemm<Tcomp>(cublas_handle, 'n', 'n', this->nactu(), this->nactu(),
                    this->nactu(), 1.0f, d_tmp->get_data(), this->nactu(),
                    d_M1->get_data(), this->nactu(), 0.0f, d_tmp2->get_data(),
                    this->nactu());

  // 4. SVdec(C',A)
  carma_magma_syevd<Tcomp>('V', d_tmp2, h_eigenvals);

  // 5. M = U * D⁻¹
  carma_dgmm<Tcomp>(cublas_handle, CUBLAS_SIDE_RIGHT, this->nactu(),
                    this->nactu(), d_geocov->get_data(), this->nactu(),
                    d_eigenvals_inv->get_data(), 1, d_tmp->get_data(),
                    this->nactu());

  // 6. B = M * A;
  carma_gemm<Tcomp>(cublas_handle, 'n', 'n', this->nactu(), this->nactu(),
                    this->nactu(), 1.0f, d_tmp->get_data(), this->nactu(),
                    d_tmp2->get_data(), this->nactu(), 0.0f,
                    d_KLbasis->get_data(), this->nactu());

  delete d_M1;
  delete d_tmp;
  delete d_tmp2;
  delete d_eigenvals;
  delete d_eigenvals_sqrt;
  delete d_eigenvals_inv;
  delete h_eigenvals;
  delete h_eigenvals_sqrt;
  delete h_eigenvals_inv;

  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int sutra_controller_mv<Tcomp, Tout>::load_covmat(Tcomp *covmat) {
  this->current_context->set_active_device(this->device, 1);
  this->d_covmat->host2device(covmat);
  return EXIT_SUCCESS;
}
template <typename Tcomp, typename Tout>
int sutra_controller_mv<Tcomp, Tout>::load_klbasis(Tcomp *klbasis) {
  this->current_context->set_active_device(this->device, 1);
  this->d_KLbasis->host2device(klbasis);
  return EXIT_SUCCESS;
}

// Florian features
template <typename Tcomp, typename Tout>
int sutra_controller_mv<Tcomp, Tout>::build_cmat(Tcomp cond) {
  this->current_context->set_active_device(this->device, 1);
  if (this->h_Cmmeigenvals != 0L) delete this->h_Cmmeigenvals;

  long Nactu = this->d_Cphim->get_dims()[1];

  // (Cmm + Cn)⁻¹
  add_md(this->d_Cmm->get_data(), this->d_Cmm->get_data(),
         this->d_noisemat->get_data(), this->nslope(),
         this->current_context->get_device(this->device));
  // invgen(this->d_Cmm,/*(Tcomp)(this->nslope()-this->nactu())*/200.0f,0);
  long dims_data1[2] = {1, this->nslope()};
  this->h_Cmmeigenvals = new CarmaHostObj<Tcomp>(dims_data1, MA_PAGELOCK);

  invgen(this->d_Cmm, this->h_Cmmeigenvals, cond);

  // Cphim * (Cmm + Cn)⁻¹
  carma_gemm(cublas_handle, 'n', 'n', Nactu, this->nslope(), this->nslope(),
             1.0f, this->d_Cphim->get_data(), Nactu, this->d_Cmm->get_data(),
             this->nslope(), 0.0f, d_cmat->get_data(), this->nactu());

  return EXIT_SUCCESS;
}
template <typename Tcomp, typename Tout>
int sutra_controller_mv<Tcomp, Tout>::filter_cmat(Tcomp cond) {
  this->current_context->set_active_device(this->device, 1);

  if (this->h_eigenvals != 0L) delete this->h_eigenvals;
  long Nactu = this->d_Cphim->get_dims()[1];
  if (Nactu < this->nactu()) {
    long *dims_data = new long[3];
    dims_data[0] = 2;

    dims_data[1] = this->nslope();
    dims_data[2] = 2;
    CarmaObj<Tcomp> *d_M =
        new CarmaObj<Tcomp>(this->current_context, dims_data);
    dims_data[1] = Nactu;
    dims_data[2] = 2;
    CarmaObj<Tcomp> *d_TT2ho =
        new CarmaObj<Tcomp>(this->current_context, dims_data);
    dims_data[1] = 2;
    dims_data[2] = this->nslope();
    CarmaObj<Tcomp> *d_M1 =
        new CarmaObj<Tcomp>(this->current_context, dims_data);
    dims_data[2] = 2;
    CarmaObj<Tcomp> *d_tmp3 =
        new CarmaObj<Tcomp>(this->current_context, dims_data);

    // Imat decomposition TT
    dims_data[1] = Nactu;
    dims_data[2] = Nactu;
    CarmaObj<Tcomp> *d_tmp2 =
        new CarmaObj<Tcomp>(this->current_context, dims_data);

    // Dm⁻¹
    carma_gemm(cublas_handle, 't', 'n', Nactu, Nactu, this->nslope(), 1.0f,
               d_imat->get_data(), this->nslope(), d_imat->get_data(),
               this->nslope(), 0.0f, d_tmp2->get_data(), Nactu);

    long dims_data1[2] = {1, Nactu};
    this->h_eigenvals = new CarmaHostObj<Tcomp>(dims_data1, MA_PAGELOCK);
    invgen(d_tmp2, this->h_eigenvals, cond);

    dims_data[1] = Nactu;
    dims_data[2] = this->nslope();
    CarmaObj<Tcomp> *d_Dm1 =
        new CarmaObj<Tcomp>(this->current_context, dims_data);
    carma_gemm(cublas_handle, 'n', 't', Nactu, this->nslope(), Nactu, 1.0f,
               d_tmp2->get_data(), Nactu, d_imat->get_data(), this->nslope(),
               0.0f, d_Dm1->get_data(), Nactu);

    delete d_tmp2;

    // TT2ho = Dm⁻¹ * Dtt
    carma_gemm(cublas_handle, 'n', 'n', Nactu, 2, this->nslope(), 1.0f,
               d_Dm1->get_data(), Nactu,
               d_imat->get_data_at(this->nslope() * (Nactu)), this->nslope(),
               0.0f, d_TT2ho->get_data(), Nactu);

    delete d_Dm1;

    // M = Dm * TT2ho
    carma_gemm(cublas_handle, 'n', 'n', this->nslope(), 2, Nactu, 1.0f,
               d_imat->get_data(), this->nslope(), d_TT2ho->get_data(), Nactu,
               0.0f, d_M->get_data(), this->nslope());

    // M⁻¹
    carma_gemm(cublas_handle, 't', 'n', 2, 2, this->nslope(), 1.0f,
               d_M->get_data(), this->nslope(), d_M->get_data(), this->nslope(),
               0.0f, d_tmp3->get_data(), 2);
    invgen(d_tmp3, 0.0f, 0);

    carma_gemm(cublas_handle, 'n', 't', 2, this->nslope(), 2, 1.0f,
               d_tmp3->get_data(), 2, d_M->get_data(), this->nslope(), 0.0f,
               d_M1->get_data(), 2);

    // M*M⁻¹
    dims_data[1] = this->nslope();
    dims_data[2] = this->nslope();
    CarmaObj<Tcomp> *d_Ftt =
        new CarmaObj<Tcomp>(this->current_context, dims_data);
    carma_gemm(cublas_handle, 'n', 'n', this->nslope(), this->nslope(), 2, 1.0f,
               d_M->get_data(), this->nslope(), d_M1->get_data(), 2, 0.0f,
               d_Ftt->get_data(), this->nslope());

    // TTcomp filter
    TT_filt(d_Ftt->get_data(), this->nslope(),
            this->current_context->get_device(this->device));

    // cmat without TT
    dims_data[1] = Nactu;
    dims_data[2] = this->nslope();
    CarmaObj<Tcomp> *d_cmat_tt =
        new CarmaObj<Tcomp>(this->current_context, dims_data);

    carma_gemm(cublas_handle, 'n', 'n', Nactu, this->nslope(), this->nslope(),
               1.0f, d_cmat->get_data(), this->nactu(), d_Ftt->get_data(),
               this->nslope(), 0.0f, d_cmat_tt->get_data(), Nactu);

    delete d_Ftt;

    // Fill CMAT
    fill_cmat(this->d_cmat->get_data(), d_cmat_tt->get_data(), d_M1->get_data(),
              this->nactu(), this->nslope(),
              this->current_context->get_device(this->device));

    delete d_M;
    delete d_tmp3;
    delete d_cmat_tt;
    delete d_TT2ho;
    delete d_M1;
  }
  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int sutra_controller_mv<Tcomp, Tout>::build_cmat(const char *dmtype,
                                                 char *method) {
  Tcomp one = 1.;
  Tcomp zero = 0.;

  this->current_context->set_active_device(this->device, 1);
  if (strcmp(method, "inv") == 0) {
    //  R = (Dt*Cn⁻¹*D + Cphi⁻¹)⁻¹*Dt*Cn⁻¹

    CarmaObj<Tcomp> *d_tmp;
    CarmaObj<Tcomp> *d_tmp2;
    CarmaObj<Tcomp> *d_tmp3;
    long *dims_data2 = new long[3];
    dims_data2[0] = 2;
    dims_data2[1] = this->nactu();
    dims_data2[2] = this->nactu();
    d_tmp = new CarmaObj<Tcomp>(this->current_context, dims_data2);
    d_tmp2 = new CarmaObj<Tcomp>(this->current_context, dims_data2);
    dims_data2[1] = this->nslope();
    dims_data2[2] = this->nactu();
    d_tmp3 = new CarmaObj<Tcomp>(this->current_context, dims_data2);

    carma_dgmm(cublas_handle, CUBLAS_SIDE_LEFT, this->nslope(), this->nactu(),
               d_imat->get_data(), this->nslope(), d_noisemat->get_data(), 1,
               d_tmp3->get_data(), this->nslope());
    carma_gemm(cublas_handle, 't', 'n', this->nactu(), this->nactu(),
               this->nslope(), one, d_tmp3->get_data(), this->nslope(),
               d_imat->get_data(), this->nslope(), zero, d_tmp2->get_data(),
               this->nactu());
    delete d_tmp3;

    carma_geam(cublas_handle, 'n', 'n', this->nactu(), this->nactu(), one,
               d_tmp2->get_data(), this->nactu(), one, d_covmat->get_data(),
               this->nactu(), d_tmp->get_data(), this->nactu());
    delete d_tmp2;

    carma_magma_potri<Tcomp>(d_tmp);

    dims_data2[1] = this->nactu();
    dims_data2[2] = this->nslope();
    d_tmp2 = new CarmaObj<Tcomp>(this->current_context, dims_data2);
    carma_gemm(cublas_handle, 'n', 't', this->nactu(), this->nslope(),
               this->nactu(), one, d_tmp->get_data(), this->nactu(),
               d_imat->get_data(), this->nslope(), zero, d_tmp2->get_data(),
               this->nactu());
    delete d_tmp;

    carma_dgmm(cublas_handle, CUBLAS_SIDE_RIGHT, this->nactu(), this->nslope(),
               d_tmp2->get_data(), this->nactu(), d_noisemat->get_data(), 1,
               d_cmat->get_data(), this->nactu());

    delete d_tmp2;
  }

  else if (strcmp(method, "n") == 0) {
    //  R = Cphi*Dt*(D*Cphi*Dt + Cn)⁻¹

    CarmaObj<Tcomp> *d_tmp;
    CarmaObj<Tcomp> *d_tmp2;
    CarmaObj<Tcomp> *d_tmp3;
    CarmaObj<Tcomp> *d_tmp4;
    long *dims_data2 = new long[3];
    dims_data2[0] = 2;
    dims_data2[1] = this->nslope();
    dims_data2[2] = this->nslope();
    d_tmp = new CarmaObj<Tcomp>(this->current_context, dims_data2);
    d_tmp2 = new CarmaObj<Tcomp>(this->current_context, dims_data2);
    d_tmp4 = new CarmaObj<Tcomp>(this->current_context, dims_data2);
    //    CarmaObj<Tcomp> *d_U;
    //    d_U = new CarmaObj<Tcomp>(this->current_context, dims_data2);
    dims_data2[1] = this->nslope();
    dims_data2[2] = this->nactu();
    d_tmp3 = new CarmaObj<Tcomp>(this->current_context, dims_data2);
    long *dims_data = new long[2];
    dims_data[0] = 1;
    dims_data[1] = this->nactu();
    // CarmaHostObj<Tcomp> *h_eigenvals = new CarmaHostObj<Tcomp>(dims_data,
    // MA_PAGELOCK); CarmaObj<Tcomp> *d_eigenvals = new
    // CarmaObj<Tcomp>(this->current_context, dims_data);

    carma_gemm(cublas_handle, 'n', 'n', this->nslope(), this->nactu(),
               this->nactu(), one, d_imat->get_data(), this->nslope(),
               d_covmat->get_data(), this->nactu(), zero, d_tmp3->get_data(),
               this->nslope());
    carma_gemm(cublas_handle, 'n', 't', this->nslope(), this->nslope(),
               this->nactu(), one, d_tmp3->get_data(), this->nslope(),
               d_imat->get_data(), this->nslope(), zero, d_tmp2->get_data(),
               this->nslope());
    delete d_tmp3;
    add_md(d_tmp4->get_data(), d_tmp4->get_data(), d_noisemat->get_data(),
           this->nslope(), this->current_context->get_device(this->device));
    carma_geam(cublas_handle, 'n', 'n', this->nslope(), this->nslope(), one,
               d_tmp2->get_data(), this->nslope(), one, d_tmp4->get_data(),
               this->nslope(), d_tmp->get_data(), this->nslope());
    delete d_tmp2;
    delete d_tmp4;

    carma_magma_potri<Tcomp>(d_tmp);

    dims_data2[1] = this->nactu();
    dims_data2[2] = this->nslope();
    d_tmp2 = new CarmaObj<Tcomp>(this->current_context, dims_data2);
    carma_gemm(cublas_handle, 't', 'n', this->nactu(), this->nslope(),
               this->nslope(), one, d_imat->get_data(), this->nslope(),
               d_tmp->get_data(), this->nslope(), zero, d_tmp2->get_data(),
               this->nactu());
    delete d_tmp;
    carma_gemm(cublas_handle, 'n', 'n', this->nactu(), this->nslope(),
               this->nactu(), one, d_covmat->get_data(), this->nactu(),
               d_tmp2->get_data(), this->nactu(), zero, d_cmat->get_data(),
               this->nactu());

    delete d_tmp2;
    //    delete d_U;
    delete[] dims_data;
    delete[] dims_data2;
  } else {
  }  // y_error("Specify the computation method for mv : inv or n \n");}
  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int sutra_controller_mv<Tcomp, Tout>::frame_delay() {
  // here we place the content of d_centroids into cenbuf and get
  // the actual centroid frame for error computation depending on delay value

  this->current_context->set_active_device(this->device, 1);
  if ((int)this->delay > 0) {
    for (int cc = 0; cc < this->delay; cc++)
      shift_buf(&((this->d_cenbuff->get_data())[cc * this->nslope()]), 1,
                this->nslope(),
                this->current_context->get_device(this->device));

    carma_safe_call(cudaMemcpy(
        &(this->d_cenbuff->get_data()[(int)this->delay * this->nslope()]),
        this->d_centroids->get_data(), sizeof(Tcomp) * this->nslope(),
        cudaMemcpyDeviceToDevice));

    carma_safe_call(
        cudaMemcpy(this->d_centroids->get_data(), this->d_cenbuff->get_data(),
                   sizeof(Tcomp) * this->nslope(), cudaMemcpyDeviceToDevice));
  }

  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int sutra_controller_mv<Tcomp, Tout>::comp_com() {
  this->current_context->set_active_device(this->device, 1);

  // this->frame_delay();

  // POLC equations

  carma_geam<Tcomp>(
      cublas_handle, 'n', 'n', this->nactu(), 1, (Tcomp)(this->a),
      this->d_com1->get_data(), this->nactu(), this->b,
      this->d_com->get_data(), this->nactu(), *d_compbuff, this->nactu());
  carma_gemv<Tcomp>(cublas_handle, 'n', this->nslope(), this->nactu(), 1.0f,
                    *d_imat, this->nslope(), *d_compbuff, 1, 0.0f, *d_compbuff2,
                    1);
  carma_geam<Tcomp>(cublas_handle, 'n', 'n', this->nslope(), 1, 1.0f,
                    *this->d_centroids, this->nslope(), -1.0f, *d_compbuff2,
                    this->nslope(), *d_olmeas, this->nslope());

  int nstreams = this->streams->get_nb_streams();
  if (nstreams > 1) {
    Tcomp alpha = -1.0f;
    Tcomp beta = 0.0f;

    for (int i = 0; i < nstreams; i++) {
      int istart1 =
          i * this->d_cmat->get_dims(2) * this->d_cmat->get_dims(1) / nstreams;
      int istart2 = i * this->d_cmat->get_dims(1) / nstreams;

      cublasSetStream(cublas_handle, this->streams->get_stream(i));

      cublasOperation_t trans = carma_char2cublas_operation('n');

      carma_checkCublasStatus(
          cublasSgemv(cublas_handle, trans, this->d_cmat->get_dims(1) / nstreams,
                      this->d_cmat->get_dims(2), &alpha,
                      &((this->d_cmat->get_data())[istart1]),
                      this->d_cmat->get_dims(1) / nstreams, *d_olmeas, 1, &beta,
                      &((this->d_err->get_data())[istart2]), 1));
    }

    this->streams->wait_all_streams();

  } else {
    // compute error
    this->d_err->gemv('n', -1.0f, this->d_cmat, this->d_cmat->get_dims(1),
                      d_olmeas, 1, 0.0f, 1);  // POLC --> d_olmeas
  }
  /*
   mult_int(this->d_com->get_data(), this->d_err->get_data(),
   this->d_gain->get_data(), this->gain, this->nactu(),
   this->current_context->get_device(device));*/

  carma_geam<Tcomp>(cublas_handle, 'n', 'n', this->nactu(), 1, this->gain,
                    *d_err, this->nactu(), 1.0f - this->gain, *this->d_com1,
                    this->nactu(), *this->d_com, this->nactu());

  return EXIT_SUCCESS;
}

template class sutra_controller_mv<float, float>;
template class sutra_controller_mv<float, uint16_t>;
