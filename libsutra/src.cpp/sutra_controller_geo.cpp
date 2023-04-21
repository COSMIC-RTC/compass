// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_controller_geo.cpp
//! \ingroup   libsutra
//! \class     sutra_controller_geo
//! \brief     this class provides the controller_geo features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.3
//! \date      2022/01/24

#include <sutra_controller_geo.h>

template <typename T, typename Tout>
sutra_controller_geo<T, Tout>::sutra_controller_geo(CarmaContext *context,
                                                    long nactu, long Nphi,
                                                    float delay, SutraDms *dms,
                                                    int *idx_dms, int ndm, int *idx_centro, int ncentro,
                                                    bool wfs_direction)
    : SutraController<T, Tout>(context, 0, nactu, 0.0f, dms, idx_dms, ndm, idx_centro, ncentro) {
  this->gain = 0.0f;
  this->Nphi = Nphi;

  //	long dims_data2[3];
  //	dims_data2[0] = 2;
  //	dims_data2[1] = nactu;
  //	dims_data2[2] = Nphi;
  // this->d_proj = new CarmaObj<float>(this->current_context, dims_data2);
  this->d_proj = 0L;
  this->d_geocov = 0L;
  this->d_IFsparse = 0L;
  this->d_geocovTT = 0L;
  this->d_TT = 0L;
  this->d_phif = 0L;
  //	this->d_Btt = 0L;
  /*
  if (delay > 0) {
      dims_data2[1] = Nphi;
      dims_data2[2] = delay + 1;
      this->d_cenbuff = new CarmaObj<float>(this->current_context, dims_data2);
    }
    */
  this->Ntt = 0;

  long dims_data1[2];
  dims_data1[0] = 1;
  dims_data1[1] = nactu;
  this->d_gain = new CarmaObj<T>(this->current_context, dims_data1);
  this->d_compfloat = new CarmaObj<float>(this->current_context, dims_data1);
  this->d_compdouble = new CarmaObj<double>(this->current_context, dims_data1);
  dims_data1[1] = Nphi;
  this->d_phi = new CarmaObj<double>(this->current_context, dims_data1);
  this->d_indx_pup = new CarmaObj<int>(this->current_context, dims_data1);
  if (wfs_direction)
    this->d_indx_mpup = new CarmaObj<int>(this->current_context, dims_data1);
  else
    this->d_indx_mpup = 0L;
}

template <typename T, typename Tout>
sutra_controller_geo<T, Tout>::~sutra_controller_geo() {
  this->current_context->set_active_device(this->device, 1);
  delete this->d_proj;
  delete this->d_gain;
  delete this->d_indx_pup;
  delete this->d_phi;
  delete this->d_compfloat;
  delete this->d_compdouble;
  if (this->Ntt) {
    delete this->d_TT;
    delete this->d_geocovTT;
    delete this->d_phif;
  }
}

template <typename T, typename Tout>
string sutra_controller_geo<T, Tout>::get_type() {
  return "geo";
}

template <typename T, typename Tout>
int sutra_controller_geo<T, Tout>::load_mgain(T *mgain) {
  this->d_gain->host2device(mgain);
  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int sutra_controller_geo<T, Tout>::load_Btt(T *Btt_pzt, T *Btt_TT) {
  // the Btt given is Btt*Btt.transpose because of computation needs
  /*
  long dims_data[3] = {2,n,m};
  if(this->d_geocov != 0L)
        delete this->d_geocov;
  this->d_geocov = new CarmaObj<T>(this->current_context, dims_data);
  */
  this->d_geocov->host2device(Btt_pzt);
  this->d_geocovTT->host2device(Btt_TT);

  return EXIT_SUCCESS;
}
template <typename T, typename Tout>
int sutra_controller_geo<T, Tout>::init_proj(SutraDms *dms, int *indx_dm,
                                             T *unitpervolt, int *indx_pup) {
  this->current_context->set_active_device(this->device, 1);
  long dims_data[3] = {2, this->Nphi, this->nactu()};
  CarmaObj<T> d_IF(this->current_context, dims_data);
  dims_data[1] = this->nactu();
  CarmaObj<T> d_tmp(this->current_context, dims_data);
  long tmp_dim = this->Nphi * this->d_dmseen.size();
  long dims_data1[2] = {1, tmp_dim};
  CarmaObj<int> d_indx(this->current_context, dims_data1, indx_dm);

  this->d_indx_pup->host2device(indx_pup);

  // Get influence functions in d_IF
  int indx_start = 0;
  int ind = 0;
  vector<SutraDm *>::iterator p;
  p = this->d_dmseen.begin();
  while (p != this->d_dmseen.end()) {
    SutraDm *dm = *p;
    dm->get_IF<T>(d_IF.get_data_at(indx_start * this->Nphi),
                  d_indx.get_data_at(this->Nphi * ind), this->Nphi,
                  1.0f /*unitpervolt[ind]*/);
    indx_start += dm->nactus;
    ind++;
    p++;
  }

  // d_tmp = (transpose(d_IF) * d_IF)⁻¹
  carma_gemm<T>(this->cublas_handle(), 't', 'n', this->nactu(), this->nactu(),
                this->Nphi, 1.0f, d_IF.get_data(), d_IF.get_dims()[1],
                d_IF.get_data(), d_IF.get_dims()[1], 0.0f, d_tmp.get_data(),
                d_tmp.get_dims()[1]);
  carma_potr_inv(&d_tmp);
  // invgen(d_tmp,1000.0f,1);

  // d_proj = d_tmp * transpose(d_IF)
  carma_gemm<T>(this->cublas_handle(), 'n', 't', this->nactu(), this->Nphi,
                this->nactu(), 1.0f, d_tmp.get_data(), d_tmp.get_dims()[1],
                d_IF.get_data(), d_IF.get_dims()[1], 0.0f,
                this->d_proj->get_data(), this->d_proj->get_dims()[1]);

  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int sutra_controller_geo<T, Tout>::init_proj_sparse(
    SutraDms *dms, int *indx_dm, T *unitpervolt, int *indx_pup, int *indx_mpup,
    bool roket) {
  this->current_context->set_active_device(this->device, 1);
  vector<SutraDm *>::iterator p;
  if (roket) {
    this->Ntt = 0;
    p = this->d_dmseen.begin();
    while (p != this->d_dmseen.end()) {
      SutraDm *dm = *p;
      if (dm->type == "tt") this->Ntt += 1;
      p++;
    }
  }

  int Npzt = this->d_dmseen.size() - this->Ntt;
  CarmaSparseObj<double> *d_IFi[Npzt];
  long tmp_dim = this->Nphi * this->d_dmseen.size();
  long dims_data1[2] = {1, tmp_dim};
  CarmaObj<int> d_indx(this->current_context, dims_data1, indx_dm);

  this->d_indx_pup->host2device(indx_pup);
  if (this->d_indx_mpup != 0L) this->d_indx_mpup->host2device(indx_mpup);

  // Get influence functions of the DM #ind in d_IFi
  int indx_start = 0;
  int ind = 0;
  int nnz = 0;
  int NNZ[Npzt];
  int Nact[Npzt];

  p = this->d_dmseen.begin();
  while (p != this->d_dmseen.end() - this->Ntt) {
    SutraDm *dm = *p;
    dm->get_IF_sparse<double>(d_IFi[ind], d_indx.get_data_at(this->Nphi * ind),
                              this->Nphi, 1.0f, 1);
    dm->reset_shape();
    NNZ[ind] = d_IFi[ind]->nz_elem;
    Nact[ind] = dm->nactus;
    nnz += d_IFi[ind]->nz_elem;
    indx_start += dm->nactus;
    ind++;
    p++;
  }
  // Create global d_IF_sparse from array of d_IFi
  long dims_data[2] = {1, nnz};
  CarmaObj<double> d_val(this->current_context, dims_data);
  CarmaObj<int> d_col(this->current_context, dims_data);
  dims_data[1] = (this->nactu() - 2 * this->Ntt) + 1;
  CarmaObj<int> d_row(this->current_context, dims_data);
  int cpt[Npzt];
  int nact = 0;
  cpt[0] = 0;
  p = this->d_dmseen.begin();

  for (int i = 0; i < Npzt; i++) {
    SutraDm *dm = *p;
    carma_safe_call(cudaMemcpyAsync(d_val.get_data_at(cpt[i]), d_IFi[i]->d_data,
                                  sizeof(double) * d_IFi[i]->nz_elem,
                                  cudaMemcpyDeviceToDevice));
    carma_safe_call(cudaMemcpyAsync(d_col.get_data_at(cpt[i]), d_IFi[i]->d_colind,
                                  sizeof(int) * d_IFi[i]->nz_elem,
                                  cudaMemcpyDeviceToDevice));
    if (i == 0)
      carma_safe_call(cudaMemcpyAsync(d_row.get_data(), d_IFi[i]->d_rowind,
                                    sizeof(int) * (dm->nactus + 1),
                                    cudaMemcpyDeviceToDevice));
    else
      carma_safe_call(cudaMemcpyAsync(
          d_row.get_data_at(nact + 1), &(d_IFi[i]->d_rowind[1]),
          sizeof(int) * (dm->nactus), cudaMemcpyDeviceToDevice));
    cpt[i + 1] = cpt[i] + d_IFi[i]->nz_elem;
    nact += dm->nactus;
    p++;
    delete d_IFi[i];
  }

  if (Npzt > 1) {
    dims_data[1] = Npzt;
    CarmaObj<int> d_NNZ(this->current_context, dims_data);
    CarmaObj<int> d_nact(this->current_context, dims_data);
    d_NNZ.host2device(NNZ);
    d_nact.host2device(Nact);

    adjust_csr_index(d_row.get_data_at(1), d_NNZ.get_data(), d_nact.get_data(),
                     this->nactu() - 2 * this->Ntt, Nact[0],
                     this->current_context->get_device(this->device));
  }

  long dims_data2[3] = {2, (this->nactu() - 2 * this->Ntt), this->Nphi};
  this->d_IFsparse = new CarmaSparseObj<double>(
      this->current_context, dims_data2, d_val.get_data(), d_col.get_data(),
      d_row.get_data(), nnz, false);

  // d_geocov = (transpose(d_IF) * d_IF)⁻¹
  CarmaSparseObj<double> *d_tmp =
      new CarmaSparseObj<double>(this->current_context);
  dims_data2[2] = (this->nactu() - 2 * this->Ntt);
  this->d_geocov = new CarmaObj<T>(this->current_context, dims_data2);
  CarmaObj<double> *d_tmp2 =
      new CarmaObj<double>(this->current_context, this->d_geocov->get_dims());

  carma_gemm<double>(cusparse_handle(), 'n', 't', this->d_IFsparse,
                     this->d_IFsparse, d_tmp);
  carma_csr2dense<double>(d_tmp, d_tmp2->get_data());
  double_to_float(d_tmp2->get_data(), this->d_geocov->get_data(),
                this->d_geocov->get_nb_elements(),
                this->current_context->get_device(this->device));

  this->d_geocov->scale(1.0f / this->Nphi, 1);
  carma_potr_inv(d_geocov);
  // invgen(d_geocov,2.0f,0);
  delete d_tmp;
  delete d_tmp2;
  if (this->Ntt) {
    dims_data2[1] = this->Nphi;     // 2*this->Ntt;
    dims_data2[2] = 2 * this->Ntt;  // this->Nphi;
    this->d_TT = new CarmaObj<T>(this->current_context, dims_data2);
    dims_data2[1] = 2 * this->Ntt;
    this->d_geocovTT = new CarmaObj<T>(this->current_context, dims_data2);
    dims_data1[1] = this->Nphi;
    this->d_phif = new CarmaObj<T>(this->current_context, dims_data1);

    p = this->d_dmseen.begin();
    ind = 0;
    int ind2 = 0;
    while (p != this->d_dmseen.end()) {
      SutraDm *dm = *p;
      if (dm->type == "tt") {
        // dm->get_IF(this->d_TT->get_data_at(ind*Nphi),
        // d_indx.get_data_at(this->Nphi * ind2), this->Nphi, 1.0f);
        for (int i = 0; i < dm->nactus; i++) {
          dm->comp_oneactu(i, 1.0f);

          getIF<T>(this->d_TT->get_data_at(ind * this->Nphi),
                   dm->d_shape->d_screen->get_data(),
                   d_indx.get_data_at(ind2 * this->Nphi), this->Nphi, 0,
                   dm->nactus, 1,
                   this->current_context->get_device(this->device));
          dm->reset_shape();

          ind++;
        }
      }
      ind2++;
      p++;
    }

    carma_gemm(this->cublas_handle(), 't', 'n', 2 * this->Ntt, 2 * this->Ntt,
               this->Nphi, 1.0f / this->Nphi, this->d_TT->get_data(), this->Nphi,
               this->d_TT->get_data(), this->Nphi, 0.0f,
               this->d_geocovTT->get_data(), 2 * this->Ntt);

    T *tmp;
    tmp = (T *)malloc(this->d_geocovTT->get_nb_elements() * sizeof(T));
    this->d_geocovTT->device2host(tmp);
    tmp[0] = 1.0f / tmp[0];
    tmp[3] = 1.0f / tmp[3];
    tmp[1] = 0.0f;
    tmp[2] = 0.0f;
    this->d_geocovTT->host2device(tmp);
    delete tmp;
  }

  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int sutra_controller_geo<T, Tout>::comp_dphi(SutraSource *target,
                                             bool wfs_direction) {
  this->current_context->set_active_device(this->device, 1);
  // Get the target phase in the pupil
  if (wfs_direction && this->d_indx_mpup == 0L) {
    DEBUG_TRACE("controller geo has not been initialized for wfs direction");
    return EXIT_FAILURE;
  }
  if (wfs_direction)
    get_pupphase(this->d_phi->get_data(), target->d_phase->d_screen->get_data(),
                 this->d_indx_mpup->get_data(), this->Nphi,
                 this->current_context->get_device(this->device));
  else
    get_pupphase(this->d_phi->get_data(), target->d_phase->d_screen->get_data(),
                 this->d_indx_pup->get_data(), this->Nphi,
                 this->current_context->get_device(this->device));

  remove_avg(this->d_phi->get_data(), this->d_phi->get_nb_elements(),
             this->current_context->get_device(this->device));

  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int sutra_controller_geo<T, Tout>::comp_com() {
  // Project the phase on the actuators
  /*
  //Dense version
  carma_gemv(this->cublas_handle(),'n', this->nactu(), this->Nphi, -1.0f,
  this->d_proj->get_data(),this->d_proj->get_dims()[1], this->d_phi->get_data(),1,
  0.0f, this->d_com->get_data(),1);
  */
  this->current_context->set_active_device(this->device, 1);

  // Sparse version
  carma_gemv(cusparse_handle(), 'n', 1.0 / this->Nphi, this->d_IFsparse,
             this->d_phi->get_data(), 0.0, this->d_compdouble->get_data());
  double_to_float(this->d_compdouble->get_data(), this->d_compfloat->get_data(),
                this->nactu(), this->current_context->get_device(this->device));
  // If we are in error budget case, d_geocov is Btt*Btt.transpose
  carma_gemv(this->cublas_handle(), 'n', this->nactu() - 2 * this->Ntt,
             this->nactu() - 2 * this->Ntt, -1.0f, this->d_geocov->get_data(),
             this->d_geocov->get_dims()[1], this->d_compfloat->get_data(), 1,
             0.0f, this->d_com->get_data(), 1);
  if (this->Ntt) {
    double_to_float(this->d_phi->get_data(), this->d_phif->get_data(), this->Nphi,
                  this->current_context->get_device(this->device));
    carma_gemv(this->cublas_handle(), 't', this->d_TT->get_dims(1),
               this->d_TT->get_dims(2), 1.0f / this->Nphi, this->d_TT->get_data(),
               this->d_TT->get_dims(1), this->d_phif->get_data(), 1, 0.0f,
               this->d_compfloat->get_data(), 1);
    carma_gemv(this->cublas_handle(), 'n', 2 * this->Ntt, 2 * this->Ntt, -1.0f,
               this->d_geocovTT->get_data(), 2 * this->Ntt,
               this->d_compfloat->get_data(), 1, 0.0f,
               this->d_com->get_data_at(this->d_IFsparse->get_dims(1)), 1);
  }

  return EXIT_SUCCESS;
}

template class sutra_controller_geo<float, float>;
template class sutra_controller_geo<float, uint16_t>;
