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

//! \file      sutra_dm.cpp
//! \ingroup   libsutra
//! \class     SutraDm
//! \brief     this class provides the dm features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#include <stdio.h>
#include <sutra_dm.hpp>

//#include <sutra_dm.cuh>

SutraDms::SutraDms() {}

SutraDms::~SutraDms() {
  for (std::vector<SutraDm *>::iterator it = this->d_dms.begin();
       this->d_dms.end() != it; it++) {
    delete *it;
  }
  this->d_dms.clear();
}

int32_t SutraDms::add_dm(CarmaContext *context, const char *type, float alt,
                      int64_t dim, int64_t nactus, int64_t influsize, int64_t ninflupos,
                      int64_t n_npoints, float push4imat, int64_t nord, int32_t device) {
  this->insert_dm(context, type, alt, dim, nactus, influsize, ninflupos,
                  n_npoints, push4imat, nord, 0.f, 0.f, 0.f, 1.f, device, this->d_dms.size());

  return EXIT_SUCCESS;
}

int32_t SutraDms::add_dm(CarmaContext *context, const char *type, float alt,
                      int64_t dim, int64_t nactus, int64_t influsize, int64_t ninflupos,
                      int64_t n_npoints, float push4imat, int64_t nord, float dx, float dy, float thetaML, float G, int32_t device) {
  this->insert_dm(context, type, alt, dim, nactus, influsize, ninflupos,
                  n_npoints, push4imat, nord, dx, dy, thetaML, G, device, this->d_dms.size());

  return EXIT_SUCCESS;
}

int32_t SutraDms::insert_dm(CarmaContext *context, const char *type, float alt,
                         int64_t dim, int64_t nactus, int64_t influsize, int64_t ninflupos,
                         int64_t n_npoints, float push4imat, int64_t nord, float dx, float dy, float thetaML, float G, int32_t device,
                         int32_t idx) {
  d_dms.insert(d_dms.begin() + idx,
               new SutraDm(context, type, alt, dim, nactus, influsize,
                            ninflupos, n_npoints, push4imat, nord, dx, dy, thetaML, G, device));
  return EXIT_SUCCESS;
}

int32_t SutraDms::remove_dm(int32_t idx) {
  if (idx < this->d_dms.size()) {
    delete d_dms[idx];
    d_dms.erase(d_dms.begin() + idx);
  }

  else
    DEBUG_TRACE("Index exceed vector size");

  return EXIT_SUCCESS;
}

int32_t SutraDms::nact_total() {
  int32_t nact = 0;
  for (size_t idx = 0; idx < this->d_dms.size(); idx++) {
    nact += d_dms[idx]->nactus;
  }
  return nact;
}

SutraDm::SutraDm(CarmaContext *context, const char *type, float alt,
                   int64_t dim, int64_t nactus, int64_t influsize, int64_t ninflupos,
                   int64_t n_npoints, float push4imat, int64_t nord, float dx, float dy, float thetaML, float G, int32_t device) {
  this->d_influ = NULL;
  this->d_influpos = NULL;
  this->d_npoints = NULL;
  this->d_istart = NULL;
  this->d_xoff = NULL;
  this->d_yoff = NULL;

  this->d_KLbasis = NULL;
  // this->d_IFsparse = NULL;
  // this->d_comdouble = NULL;
  // this->d_shapedouble = NULL;

  this->current_context = context;
  this->device = device;
  current_context->set_active_device(device, 1);
  this->nactus = nactus;
  this->dim = dim;
  this->influsize = influsize;
  this->d_shape = new SutraPhase(context, dim);
  this->type = type;
  this->altitude = alt;
  this->push4imat = push4imat;
  this->volt_min = -1.0f;
  this->volt_max = 1.0f;
  this->val_max = uint16_t(65535);
  this->dx = dx;
  this->dy = dy;
  this->thetaML = thetaML;
  this->G = G;

  int64_t dims_data1[2];
  dims_data1[0] = 1;
  dims_data1[1] = nactus;
  this->d_com = new CarmaObj<float>(context, dims_data1);

  if (strcmp(type, "kl") != 0) {
    int64_t *dims_data3 = new int64_t[4];
    dims_data3[0] = 3;
    dims_data3[1] = influsize;
    dims_data3[2] = influsize;
    dims_data3[3] = nactus;
    this->d_influ = new CarmaObj<float>(context, dims_data3);
    delete[] dims_data3;

    dims_data1[1] = nactus;
    this->d_xoff = new CarmaObj<int32_t>(context, dims_data1);
    this->d_yoff = new CarmaObj<int32_t>(context, dims_data1);
  }

  if (strcmp(type, "pzt") == 0) {
    dims_data1[1] = ninflupos;
    this->d_influpos = new CarmaObj<int32_t>(context, dims_data1);

    dims_data1[1] = n_npoints;  // *2;
    this->d_npoints = new CarmaObj<int32_t>(context, dims_data1);
    dims_data1[1] = n_npoints + 1;
    this->d_istart = new CarmaObj<int32_t>(context, dims_data1);
  }
  if (strcmp(type, "kl") == 0) {
    this->d_kl = new SutraKL(context, influsize, ninflupos, n_npoints, nactus,
                              nord, device);
  }
}

SutraDm::~SutraDm() {
  current_context->set_active_device(device, 1);

  delete this->d_shape;
  delete this->d_com;

  if (this->d_influ != NULL) delete this->d_influ;
  if (this->d_influpos != NULL) delete this->d_influpos;
  if (this->d_npoints != NULL) delete this->d_npoints;
  if (this->d_istart != NULL) delete this->d_istart;
  if (this->d_xoff != NULL) delete this->d_xoff;
  if (this->d_yoff != NULL) delete this->d_yoff;
  if (this->d_KLbasis != NULL) delete this->d_KLbasis;
}

int32_t SutraDm::nact() { return this->nactus; }

int32_t SutraDm::set_registration(float dx, float dy, float thetaML, float G) {
  this->dx = dx;
  this->dy = dy;
  this->thetaML = thetaML;
  this->G = G;

  return EXIT_SUCCESS;
}

int32_t SutraDm::tt_loadarrays(float *influ) {
  current_context->set_active_device(device, 1);
  this->d_influ->host2device(influ);
  return EXIT_SUCCESS;
}

int32_t SutraDm::pzt_loadarrays(float *influ, int32_t *influpos, int32_t *npoints,
                             int32_t *istart, int32_t *xoff, int32_t *yoff) {
  // current_context->set_active_device(device, 1);
  this->d_influ->host2device(influ);

  this->d_xoff->host2device(xoff);
  this->d_yoff->host2device(yoff);
  this->d_istart->host2device(istart);
  this->d_influpos->host2device(influpos);
  this->d_npoints->host2device(npoints);

  return EXIT_SUCCESS;
}

int32_t SutraDm::kl_loadarrays(float *rabas, float *azbas, int32_t *ord, float *cr,
                            float *cp) {
  current_context->set_active_device(device, 1);
  this->d_kl->d_rabas->host2device(rabas);
  this->d_kl->d_azbas->host2device(azbas);
  this->d_kl->h_ord->fill_from(ord);
  this->d_kl->d_ord->host2device(ord);
  this->d_kl->d_cr->host2device(cr);
  this->d_kl->d_cp->host2device(cp);

  return EXIT_SUCCESS;
}

int32_t SutraDm::reset_shape() {
  current_context->set_active_device(device, 1);
  this->d_shape->d_screen->reset();
  return EXIT_SUCCESS;
}

int32_t SutraDm::comp_shape(uint16_t *comvec) {
  current_context->set_active_device(device, 1);
  convertToCom(comvec, this->d_com->get_data(), this->d_com->get_nb_elements(),
               this->volt_min, this->volt_max, this->val_max,
               this->current_context->get_device(this->device));
  return this->comp_shape();
}
int32_t SutraDm::comp_shape() { return this->comp_shape(this->d_com->get_data()); }

#ifdef CHEAT_CODE
int32_t SutraDm::comp_shape(float *comvec) {
  current_context->set_active_device(device, 1);
  this->reset_shape();

  int32_t nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(current_context->get_device(device),
                         this->d_shape->d_screen->get_nb_elements(), nb_blocks,
                         nb_threads);

  if (this->type == "pzt") {
    comp_dmshape(nb_threads, nb_blocks, this->d_influ->get_data(),
                 this->d_shape->d_screen->get_data(),
                 this->d_influpos->get_data(), this->d_istart->get_data(),
                 this->d_npoints->get_data(), comvec,
                 this->influsize * this->influsize,
                 this->d_shape->d_screen->get_nb_elements());
  }

  if (this->type == "tt")
    comp_fulldmshape(nb_threads, nb_blocks, this->d_influ->get_data(),
                     this->d_shape->d_screen->get_data(), this->nactus,
                     this->influsize * this->influsize, comvec,
                     this->d_shape->d_screen->get_nb_elements());

  if (this->type == "kl") {
    int32_t xoff =
        (int32_t)((this->d_shape->d_screen->get_dims()[1] - this->d_kl->dim) / 2.0f);
    int32_t yoff = xoff;
    this->d_kl->do_combi(comvec, this->d_shape->d_screen->get_data(),
                         this->d_shape->d_screen->get_dims()[1], xoff, yoff);
  }
  return EXIT_SUCCESS;
}

#else

int32_t SutraDm::comp_shape(float *comvec) {
  current_context->set_active_device(device, 1);
  this->reset_shape();

  dim3 threads(BLOCKSIZE);
  dim3 blocks(CEIL(this->d_shape->d_screen->get_nb_elements() << 2, threads.x));
  int32_t shared = 0;

  int32_t nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(current_context->get_device(device),
                         this->d_shape->d_screen->get_nb_elements(), nb_blocks,
                         nb_threads);

  if (this->type == "pzt") {
    comp_dmshape2<float>(
        this->d_shape->d_screen->get_data(), comvec, this->d_influ->get_data(),
        this->d_istart->get_data(), this->d_npoints->get_data(),
        this->d_shape->d_screen->get_nb_elements(), threads, blocks, shared);
  }

  if (this->type == "tt")
    comp_fulldmshape(nb_threads, nb_blocks, this->d_influ->get_data(),
                     this->d_shape->d_screen->get_data(), this->nactus,
                     this->influsize * this->influsize, comvec,
                     this->d_shape->d_screen->get_nb_elements());

  if (this->type == "kl") {
    int32_t xoff =
        (int32_t)((this->d_shape->d_screen->get_dims()[1] - this->d_kl->dim) / 2.0f);
    int32_t yoff = xoff;
    this->d_kl->do_combi(comvec, this->d_shape->d_screen->get_data(),
                         this->d_shape->d_screen->get_dims()[1], xoff, yoff);
  }
  return EXIT_SUCCESS;
}

#endif

int32_t SutraDm::comp_oneactu(int32_t nactu, float ampli) {
  current_context->set_active_device(device, 1);
  this->reset_shape();
  int32_t nb_threads = 0, nb_blocks = 0;
  // get_num_blocks_and_threads(this->device,this->dim * this->dim, nb_blocks,
  // nb_threads);
  get_num_blocks_and_threads(current_context->get_device(device),
                         this->influsize * this->influsize, nb_blocks, nb_threads);
  if (this->type == "pzt")
    oneactu(nb_threads, nb_blocks, this->d_influ->get_data(),
            this->d_shape->d_screen->get_data(), nactu, ampli,
            this->d_xoff->get_data(), this->d_yoff->get_data(), this->dim,
            this->influsize, this->influsize * this->influsize);
  if (this->type == "tt")
    oneactu(nb_threads, nb_blocks, this->d_influ->get_data(),
            this->d_shape->d_screen->get_data(), nactu, ampli, this->dim,
            this->influsize, this->influsize * this->influsize);
  if (this->type == "kl") {
    int32_t xoff =
        (int32_t)((this->d_shape->d_screen->get_dims()[1] - this->d_kl->dim) / 2.0f);
    int32_t yoff = xoff;
    this->d_kl->do_compute(ampli, this->d_shape->d_screen->get_data(), nactu,
                           this->d_shape->d_screen->get_dims()[1], xoff, yoff);
  }

  return EXIT_SUCCESS;
}

template <class T>
int32_t SutraDm::get_IF(T *IF, int32_t *indx_pup, int64_t nb_pts, float ampli) {
  for (int32_t i = 0; i < this->nactus; i++) {
    this->comp_oneactu(i, ampli);
    getIF<T>(IF, this->d_shape->d_screen->get_data(), indx_pup, nb_pts, i,
             this->nactus, 1, current_context->get_device(device));
  }

  this->reset_shape();

  return EXIT_SUCCESS;
}
template int32_t SutraDm::get_IF<float>(float *IF, int32_t *indx_pup, int64_t nb_pts,
                                     float ampli);
template int32_t SutraDm::get_IF<double>(double *IF, int32_t *indx_pup, int64_t nb_pts,
                                      float ampli);

template <class T>
int32_t SutraDm::get_IF_sparse(CarmaSparseObj<T> *&d_IFsparse, int32_t *indx_pup,
                            int64_t nb_pts, float ampli, int32_t puponly) {
  current_context->set_active_device(device, 1);
  int32_t nnz_tot = 0;
  float *values[this->nactus];
  int32_t *colind[this->nactus];
  int32_t NZ[this->nactus];
  int64_t dims_data2[3] = {2, 1, nb_pts};
  CarmaObj<T> d_IF(current_context, dims_data2);
  CarmaSparseObj<T> *d_IFsparse_vec;

  std::cout << "Computing IF sparse..." << std::endl;
  for (int32_t i = 0; i < this->nactus; i++) {
    // Compute and store IF for actu i in d_IF
    this->comp_oneactu(i, ampli);
    getIF<T>(d_IF.get_data(), this->d_shape->d_screen->get_data(), indx_pup,
             nb_pts, 0, this->nactus, puponly,
             this->current_context->get_device(device));
    // CUsparse d_IF
    d_IFsparse_vec = new CarmaSparseObj<T>(&d_IF);
    // Retrieve nnz, values and colind from d_IFsparse_vec, stored on CPU
    // DEBUG_TRACE("nnz : %d \n",d_IFsparse_vec->get_nonzero_elem());

    NZ[i] = d_IFsparse_vec->get_nonzero_elem();
    values[i] = (float *)malloc(NZ[i] * sizeof(T));
    colind[i] = (int32_t *)malloc(NZ[i] * sizeof(int32_t));

    carma_safe_call(cudaMemcpyAsync(values[i], d_IFsparse_vec->get_data(),
                                  sizeof(T) * NZ[i], cudaMemcpyDeviceToHost));
    carma_safe_call(cudaMemcpyAsync(colind[i], d_IFsparse_vec->d_colind,
                                  sizeof(int32_t) * NZ[i], cudaMemcpyDeviceToHost));

    nnz_tot += NZ[i];

    delete d_IFsparse_vec;
  }
  // Reconstruction of d_data, d_colind, d_rowind for IFsparse
  int64_t dims_data[2] = {1, nnz_tot};
  CarmaObj<T> d_val(current_context, dims_data);
  CarmaObj<int32_t> d_col(current_context, dims_data);
  dims_data[1] = this->nactus + 1;
  CarmaObj<int32_t> d_row(current_context, dims_data);
  int32_t cpt[this->nactus + 1];
  cpt[0] = 0;

  for (int32_t i = 0; i < this->nactus; i++) {
    carma_safe_call(cudaMemcpyAsync(d_val.get_data_at(cpt[i]), values[i],
                                  sizeof(T) * NZ[i], cudaMemcpyHostToDevice));
    carma_safe_call(cudaMemcpyAsync(d_col.get_data_at(cpt[i]), colind[i],
                                  sizeof(int32_t) * NZ[i], cudaMemcpyHostToDevice));
    cpt[i + 1] = cpt[i] + NZ[i];
  }
  carma_safe_call(cudaMemcpyAsync(d_row.get_data(), cpt,
                                sizeof(int32_t) * (this->nactus + 1),
                                cudaMemcpyHostToDevice));
  dims_data2[1] = this->nactus;

  d_IFsparse =
      new CarmaSparseObj<T>(current_context, dims_data2, d_val.get_data(),
                              d_col.get_data(), d_row.get_data(), nnz_tot, false);

  return EXIT_SUCCESS;
}
template int32_t SutraDm::get_IF_sparse<float>(
    CarmaSparseObj<float> *&d_IFsparse, int32_t *indx_pup, int64_t nb_pts,
    float ampli, int32_t puponly);
template int32_t SutraDm::get_IF_sparse<double>(
    CarmaSparseObj<double> *&d_IFsparse, int32_t *indx_pup, int64_t nb_pts,
    float ampli, int32_t puponly);

int32_t SutraDm::compute_KLbasis(float *xpos, float *ypos, int32_t *indx, int64_t dim,
                              float norm, float ampli) {
  current_context->set_active_device(device, 1);
  int64_t dims_data[3];
  dims_data[0] = 2;
  dims_data[1] = this->nactus;
  dims_data[2] = this->nactus;
  CarmaObj<float> *d_statcov =
      new CarmaObj<float>(current_context, dims_data);
  this->d_KLbasis = new CarmaObj<float>(current_context, dims_data);
  int64_t dims_data2[2];
  dims_data2[0] = 1;
  dims_data2[1] = dim;
  CarmaObj<int32_t> *d_indx =
      new CarmaObj<int32_t>(this->current_context, dims_data2);

  // Compute the statistic matrix from actuators positions & Kolmogorov
  // statistic
  dims_data2[1] = this->nactus;
  CarmaObj<float> *d_xpos = new CarmaObj<float>(current_context, dims_data2);
  CarmaObj<float> *d_ypos = new CarmaObj<float>(current_context, dims_data2);

  d_xpos->host2device(xpos);
  d_ypos->host2device(ypos);
  dm_dostatmat(d_statcov->get_data(), this->nactus, d_xpos->get_data(),
               d_ypos->get_data(), norm, current_context->get_device(device));

  delete d_xpos;
  delete d_ypos;

  // Compute and apply piston filter
  this->piston_filt(d_statcov);

  dims_data[1] = this->nactus;
  dims_data[2] = this->nactus;
  CarmaObj<float> *d_geocov = new CarmaObj<float>(current_context, dims_data);
  d_indx->host2device(indx);

  // Sparse version for geomat
  CarmaSparseObj<double> *d_IFsparsepup;
  CarmaObj<double> *d_geodouble =
      new CarmaObj<double>(current_context, d_geocov->get_dims());
  this->get_IF_sparse<double>(d_IFsparsepup, d_indx->get_data(), dim, ampli, 1);
  this->do_geomat_from_sparse(d_geodouble->get_data(), d_IFsparsepup);
  double_to_float(d_geodouble->get_data(), d_geocov->get_data(),
                d_geocov->get_nb_elements(), current_context->get_device(device));

  delete d_geodouble;
  delete d_IFsparsepup;
  // Dense version for geomat
  /*
  dims_data[1] = dim;
  dims_data[2] = this->nactus;

  CarmaObj<float> *d_IF = new CarmaObj<float>(this->current_context,
  dims_data);
  // Get influence functions of the DM
  this->get_IF(d_IF->get_data(),d_indx->get_data(),dim,ampli);
  // Compute geometric matrix (to be CUsparsed)
  this->do_geomat(d_geocov->get_data(),d_IF->get_data(),dim);
   */

  // Double diagonalisation to obtain KL basis on actuators
  this->DDiago(d_statcov, d_geocov);

  delete d_geocov;
  delete d_indx;

  return EXIT_SUCCESS;
}

template <class T>
int32_t SutraDm::do_geomat_from_sparse(T *d_geocov,
                                  CarmaSparseObj<T> *d_IFsparse) {
  current_context->set_active_device(device, 1);
  CarmaSparseObj<T> *d_tmp = new CarmaSparseObj<T>(this->current_context);

  carma_gemm<T>(cusparse_handle(), 'n', 't', d_IFsparse, d_IFsparse, d_tmp);
  carma_csr2dense<T>(d_tmp, d_geocov);

  delete d_tmp;
  return EXIT_SUCCESS;
}
template int32_t SutraDm::do_geomat_from_sparse<float>(
    float *d_geocov, CarmaSparseObj<float> *d_IFsparse);
template int32_t SutraDm::do_geomat_from_sparse<double>(
    double *d_geocov, CarmaSparseObj<double> *d_IFsparse);

int32_t SutraDm::do_geomat(float *d_geocov, float *d_IF, int64_t n_pts) {
  current_context->set_active_device(device, 1);
  carma_gemm(this->cublas_handle(), 't', 'n', this->nactus, this->nactus, n_pts,
             1.0f, d_IF, n_pts, d_IF, n_pts, 0.0f, d_geocov, this->nactus);

  return EXIT_SUCCESS;
}

int32_t SutraDm::piston_filt(CarmaObj<float> *d_statcov) {
  current_context->set_active_device(device, 1);
  int64_t Nmod = d_statcov->get_dims()[1];
  int64_t dims_data[3];
  dims_data[0] = 2;
  dims_data[1] = Nmod;
  dims_data[2] = Nmod;
  CarmaObj<float> *d_F = new CarmaObj<float>(current_context, dims_data);
  CarmaObj<float> *d_tmp = new CarmaObj<float>(current_context, dims_data);

  int32_t N = d_statcov->get_dims()[1] * d_statcov->get_dims()[1];
  fill_filtermat(d_F->get_data(), Nmod, N, current_context->get_device(device));

  carma_gemm(this->cublas_handle(), 'n', 'n', Nmod, Nmod, Nmod, 1.0f,
             d_F->get_data(), Nmod, d_statcov->get_data(), Nmod, 0.0f,
             d_tmp->get_data(), Nmod);
  carma_gemm(this->cublas_handle(), 'n', 'n', Nmod, Nmod, Nmod, 1.0f,
             d_tmp->get_data(), Nmod, d_F->get_data(), Nmod, 0.0f,
             d_statcov->get_data(), Nmod);

  delete d_tmp;
  delete d_F;

  return EXIT_SUCCESS;
}

int32_t SutraDm::DDiago(CarmaObj<float> *d_statcov, CarmaObj<float> *d_geocov) {
  current_context->set_active_device(device, 1);
  const int64_t dims_data[3] = {2, this->nactus, this->nactus};
  CarmaObj<float> *d_M1 = new CarmaObj<float>(current_context, dims_data);
  CarmaObj<float> *d_tmp = new CarmaObj<float>(current_context, dims_data);
  CarmaObj<float> *d_tmp2 = new CarmaObj<float>(current_context, dims_data);

  const int64_t dims_data2[2] = {1, this->nactus};
  CarmaObj<float> *d_eigenvals =
      new CarmaObj<float>(current_context, dims_data2);
  CarmaObj<float> *d_eigenvals_sqrt =
      new CarmaObj<float>(current_context, dims_data2);
  CarmaObj<float> *d_eigenvals_inv =
      new CarmaObj<float>(current_context, dims_data2);
  CarmaHostObj<float> *h_eigenvals =
      new CarmaHostObj<float>(dims_data2, MA_PAGELOCK);
  CarmaHostObj<float> *h_eigenvals_inv =
      new CarmaHostObj<float>(dims_data2, MA_PAGELOCK);
  CarmaHostObj<float> *h_eigenvals_sqrt =
      new CarmaHostObj<float>(dims_data2, MA_PAGELOCK);

  // 1. SVdec(geocov,U) --> Ut * geocov * U = D������
  carma_syevd<float>(SOLVER_EIG_MODE_VECTOR, d_geocov, d_eigenvals);
  d_eigenvals->device2host(h_eigenvals->get_data());
  for (int32_t i = 0; i < this->nactus; i++) {
    if(h_eigenvals->get_data()[i] < 0.f) {
      h_eigenvals->get_data()[i] = 0;
      h_eigenvals_sqrt->get_data()[i] = 0;
      h_eigenvals_inv->get_data()[i] = 0;
    }
    else {
      h_eigenvals_sqrt->get_data()[i] =
          sqrt(h_eigenvals->get_data()[i]);  // D = sqrt(D������)
      h_eigenvals_inv->get_data()[i] =
          1. /
          sqrt(h_eigenvals->get_data()[i]);  // D��������������� = 1/sqrt(D������)
    }
  }
  d_eigenvals_sqrt->host2device(h_eigenvals_sqrt->get_data());
  d_eigenvals_inv->host2device(h_eigenvals_inv->get_data());
  d_eigenvals->host2device(h_eigenvals->get_data());

  // 2. M��������������� = sqrt(eigenvals) * Ut : here, we have
  // transpose(M���������������)
  /*
  carma_dgmm<float>(this->cublas_handle(),CUBLAS_SIDE_RIGHT,this->nactus,this->nactus,
  d_geocov->get_data(), this->nactus, d_eigenvals_inv->get_data(),1,
  d_M1->get_data(), this->nactus);*/

  carma_dgmm<float>(this->cublas_handle(), CUBLAS_SIDE_RIGHT, this->nactus,
                    this->nactus, d_geocov->get_data(), this->nactus,
                    d_eigenvals_sqrt->get_data(), 1, d_M1->get_data(),
                    this->nactus);

  // 3. C' = M��������������� * statcov * M���������������t
  carma_gemm<float>(this->cublas_handle(), 't', 'n', nactus, nactus, nactus,
                    1.0f, d_M1->get_data(), nactus, d_statcov->get_data(), nactus,
                    0.0f, d_tmp->get_data(), nactus);

  carma_gemm<float>(this->cublas_handle(), 'n', 'n', nactus, nactus, nactus,
                    1.0f, d_tmp->get_data(), nactus, d_M1->get_data(), nactus,
                    0.0f, d_tmp2->get_data(), nactus);

  // 4. SVdec(C',A)
  carma_syevd<float>(SOLVER_EIG_MODE_VECTOR, d_tmp2, d_eigenvals);

  // 5. M = U * D���������������
  carma_dgmm<float>(this->cublas_handle(), CUBLAS_SIDE_RIGHT, this->nactus,
                    this->nactus, d_geocov->get_data(), this->nactus,
                    d_eigenvals_inv->get_data(), 1, d_tmp->get_data(),
                    this->nactus);

  // 6. B = M * A;
  carma_gemm<float>(this->cublas_handle(), 'n', 'n', nactus, nactus, nactus,
                    1.0f, d_tmp->get_data(), nactus, d_tmp2->get_data(), nactus,
                    0.0f, d_KLbasis->get_data(), nactus);

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
