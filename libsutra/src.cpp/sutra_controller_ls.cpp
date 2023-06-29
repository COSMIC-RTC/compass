// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_controller_ls.cpp
//! \ingroup   libsutra
//! \class     sutra_controller_ls
//! \brief     this class provides the controller_ls features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.4
//! \date      2022/01/24

#include <carma_magma.h>
#include <sutra_controller_ls.h>
#include <string>

template <typename T, typename Tout>
sutra_controller_ls<T, Tout>::sutra_controller_ls(CarmaContext *context,
                                                  long nslope, long nactu, float delay,
                                                  SutraDms *dms, int *idx_dms,
                                                  int ndm, int *idx_centro, int ncentro)
    : SutraController<T, Tout>(context, nslope, nactu, delay, dms,
                                idx_dms, ndm, idx_centro, ncentro) {
  this->d_imat = 0L;
  this->d_cmat = 0L;
  this->d_eigenvals = 0L;
  this->h_eigenvals = 0L;
  this->d_cenbuff = 0L;

  this->gain = 0.0f;

  this->is_modopti = 0;

  long dims_data1[2] = {1, 0};
  long dims_data2[3] = {2, 0, 0};

  dims_data2[1] = nslope;
  dims_data2[2] = nactu;
  this->d_imat = new CarmaObj<T>(context, dims_data2);

  dims_data2[1] = nactu;
  dims_data2[2] = nslope;
  this->d_cmat = new CarmaObj<T>(context, dims_data2);

  dims_data2[1] = dims_data2[2] = nactu;
  d_U = new CarmaObj<T>(context, dims_data2);

  dims_data1[1] = nslope < nactu ? nslope : nactu;
  this->d_eigenvals = new CarmaObj<T>(context, dims_data1);
  this->h_eigenvals = new CarmaHostObj<T>(dims_data1, MA_PAGELOCK);

  if (delay > 0) {
    dims_data2[1] = nslope;
    dims_data2[2] = (int)delay + 1;
    this->d_cenbuff = new CarmaObj<T>(context, dims_data2);
  } else
    this->d_cenbuff = 0L;
  this->d_M2V = 0L;
  this->d_S2M = 0L;
  this->d_slpol = 0L;
  this->d_Hcor = 0L;
  this->d_compbuff = 0L;
  this->d_compbuff2 = 0L;
  this->cpt_rec = 0;
  this->nrec = 0;
  this->nmodes = 0;
  this->gmin = 0.0f;
  this->gmax = 0.0f;
  this->ngain = 0;
  this->Fs = 0.0f;

  dims_data1[1] = nactu;
  this->d_err = new CarmaObj<T>(context, dims_data1);
  this->d_gain = new CarmaObj<T>(context, dims_data1);
}

template <typename T, typename Tout>
sutra_controller_ls<T, Tout>::~sutra_controller_ls() {
  this->current_context->set_active_device(this->device, 1);

  delete this->d_imat;
  delete this->d_cmat;
  delete this->d_U;
  delete this->d_eigenvals;
  delete this->h_eigenvals;

  if (this->d_cenbuff) delete this->d_cenbuff;
  delete this->d_err;
  delete this->d_gain;

  if (this->is_modopti) {
    delete this->d_M2V;
    delete this->d_S2M;
    delete this->d_slpol;
    delete this->d_Hcor;
    delete this->d_compbuff;
    delete this->d_compbuff2;
  }
}

template <typename T, typename Tout>
string sutra_controller_ls<T, Tout>::get_type() {
  return "ls";
}

template <typename T, typename Tout>
int sutra_controller_ls<T, Tout>::svdec_imat() {
  // doing U = Dt.D where D is i_mat
  T one = 1., zero = 0.;
  int nCols = this->d_imat->get_dims(2);

  this->current_context->set_active_device(this->device, 1);
  if (d_U->get_dims(1) !=
      nCols) {  // d_imat shape was changed during modal basis. Adapt yourself.
    delete d_U;
    this->d_U = new CarmaObj<T>(this->current_context,
                                 std::vector<long>{2, nCols, nCols}.data());
    delete d_eigenvals;
    delete h_eigenvals;
    this->d_eigenvals = new CarmaObj<T>(this->current_context,
                                         std::vector<long>{1, nCols}.data());
    this->h_eigenvals =
        new CarmaHostObj<T>(std::vector<long>{1, nCols}.data(), MA_PAGELOCK);
  }

  if (carma_syrk<T>(this->cublas_handle(), CUBLAS_FILL_MODE_LOWER, 't', nCols,
                    this->nslope(), one, *d_imat, this->nslope(), zero, *d_U,
                    nCols)) {
    return EXIT_FAILURE;
  }

  // we can skip this step syevd use only the lower part
  fill_sym_matrix('U', d_U->get_data(), nCols, nCols * nCols,
                  this->current_context->get_device(this->device));
  if (!carma_magma_disabled()) {

    // doing evd of U inplace
    if (carma_magma_syevd<T>(SOLVER_EIG_MODE_VECTOR, d_U, h_eigenvals) == EXIT_FAILURE) {
      // Case where MAGMA is not feeling good :-/
      return EXIT_FAILURE;
    }
    d_eigenvals->host2device(h_eigenvals->get_data());
  } else {  // CUsolver case
    // doing evd of U inplace
    if (carma_syevd<T>(SOLVER_EIG_MODE_VECTOR, d_U, d_eigenvals) == EXIT_FAILURE) {
      return EXIT_FAILURE;
    }
    d_eigenvals->device2host(h_eigenvals->get_data());
  }
  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int sutra_controller_ls<T, Tout>::set_cmat(T *cmat) {
  this->current_context->set_active_device(this->device, 1);
  this->d_cmat->host2device(cmat);
  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int sutra_controller_ls<T, Tout>::set_imat(T *imat) {
  this->current_context->set_active_device(this->device, 1);
  this->d_imat->host2device(imat);
  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int sutra_controller_ls<T, Tout>::set_modal_gains(T *mgain) {
  this->current_context->set_active_device(this->device, 1);
  this->d_gain->host2device(mgain);
  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int sutra_controller_ls<T, Tout>::build_cmat(int nfilt, bool filt_tt) {
  this->current_context->set_active_device(this->device, 1);

  long dims_data1[2] = {1, 0};
  long dims_data2[3] = {2, 0, 0};

  dims_data2[1] = dims_data2[2] = this->nactu();
  CarmaObj<T> d_tmp(this->current_context, dims_data2),
      d_tmp2(this->current_context, dims_data2);

  dims_data1[1] = this->nactu();
  CarmaObj<T> d_eigenvals_inv(this->current_context, dims_data1);
  CarmaHostObj<T> h_eigenvals_inv(dims_data1, MA_PAGELOCK);

  T one = 1., zero = 0.;

  int nb_elem = this->h_eigenvals->get_nb_elements();
  memset(h_eigenvals_inv.get_data(), 0, sizeof(T) * nb_elem);

  // filtering modes
  /*
  if (filt_tt)
    nb_elem -= 2;
  */
  for (int cc = nfilt; cc < nb_elem; cc++) {
    T eigenval = (*this->h_eigenvals)[cc];
    h_eigenvals_inv[cc] =  // 1.0f / eigenval;
        (fabs(eigenval) > 1.e-9) ? 1.0f / eigenval : 0.f;
  }
  d_eigenvals_inv.host2device(h_eigenvals_inv.get_data());

  carma_dgmm(this->cublas_handle(), CUBLAS_SIDE_RIGHT, this->nactu(),
             this->nactu(), d_U->get_data(), this->nactu(),
             d_eigenvals_inv.get_data(), one, d_tmp.get_data(), this->nactu());
  carma_gemm(this->cublas_handle(), 'n', 't', this->nactu(), this->nactu(),
             this->nactu(), one, d_tmp.get_data(), this->nactu(), d_U->get_data(),
             this->nactu(), zero, d_tmp2.get_data(), this->nactu());
  carma_gemm(this->cublas_handle(), 'n', 't', this->nactu(), this->nslope(),
             this->nactu(), one, d_tmp2.get_data(), this->nactu(),
             d_imat->get_data(), this->nslope(), zero, d_cmat->get_data(),
             this->nactu());

  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int sutra_controller_ls<T, Tout>::build_cmat(int nfilt) {
  this->current_context->set_active_device(this->device, 1);
  return this->build_cmat(nfilt, false);
}

template <typename T, typename Tout>
int sutra_controller_ls<T, Tout>::frame_delay() {
  // here we place the content of d_centroids into cenbuf and get
  // the actual centroid frame for error computation depending on delay value

  this->current_context->set_active_device(this->device, 1);
  if (this->delay > 0) {
    for (int cc = 0; cc < this->delay; cc++)
      shift_buf(this->d_cenbuff->get_data_at(cc * this->nslope()), 1,
                this->nslope(),
                this->current_context->get_device(this->device));

    carma_safe_call(cudaMemcpy(
        this->d_cenbuff->get_data_at((int)this->delay * this->nslope()),
        this->d_centroids->get_data(), sizeof(T) * this->nslope(),
        cudaMemcpyDeviceToDevice));

    carma_safe_call(
        cudaMemcpy(this->d_centroids->get_data(), this->d_cenbuff->get_data(),
                   sizeof(T) * this->nslope(), cudaMemcpyDeviceToDevice));
  }

  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int sutra_controller_ls<T, Tout>::comp_com() {
  this->current_context->set_active_device(this->device, 1);

  // this->frame_delay();

  // Modal Control Optimization
  if (this->is_modopti) {
    // Refresh when enough slopes have been recorded
    if (this->cpt_rec >= this->nrec + this->delay) {
      std::cout << "Refreshing modal gains..." << std::endl;
      modalControlOptimization();
      this->cpt_rec = 0;
    }
  }

  // INTEGRATOR
  //    T *cmat=(T*)malloc(this->d_cmat->get_nb_elements()*sizeof(T));
  //    d_cmat->device2host(cmat);
  //    DEBUG_TRACE("here %f %f %d", cmat[0], this->gain, this->open_loop);
  // compute error
  this->d_err->gemv('n', -1.0f, this->d_cmat, this->d_cmat->get_dims(1),
                    this->d_centroids, 1, 0.0f, 1);

  // apply modal gain & loop gain
  if (this->is_modopti)
    mult_int(this->d_com->get_data(), this->d_err->get_data(), this->gain,
              this->nactu(), this->current_context->get_device(this->device));
  else
    mult_int(this->d_com->get_data(), this->d_err->get_data(),
              this->d_gain->get_data(), this->gain, this->nactu(),
              this->current_context->get_device(this->device));


  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int sutra_controller_ls<T, Tout>::build_cmat_modopti() {
  this->current_context->set_active_device(this->device, 1);
  long dims_data2[3] = {2, this->nactu(), this->nmodes};
  CarmaObj<T> d_tmp(this->current_context, dims_data2);

  // Compute cmat as M2V*(modal gains)*S2M
  carma_dgmm(this->cublas_handle(), CUBLAS_SIDE_RIGHT, this->nactu(),
             this->nmodes, this->d_M2V->get_data(), this->nactu(),
             this->d_gain->get_data(), 1, d_tmp.get_data(), this->nactu());
  carma_gemm(this->cublas_handle(), 'n', 'n', this->nactu(), this->nslope(),
             this->nmodes, 1.0f, d_tmp.get_data(), this->nactu(),
             this->d_S2M->get_data(), this->nmodes, 0.0f,
             this->d_cmat->get_data(), this->nactu());

  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int sutra_controller_ls<T, Tout>::init_modalOpti(int nmodes, int nrec, T *M2V,
                                                 T gmin, T gmax, int ngain,
                                                 T Fs) {
  this->current_context->set_active_device(this->device, 1);
  this->is_modopti = 1;
  this->cpt_rec = 0;
  this->nrec = nrec;
  this->nmodes = nmodes;
  this->gmin = gmin;
  this->gmax = gmax;
  this->ngain = ngain;
  this->gain = 1.0f;
  this->Fs = Fs;
  long dims_data1[2] = {1, nmodes};
  this->d_gain = new CarmaObj<T>(this->current_context, dims_data1);
  dims_data1[1] = this->nactu();
  this->d_compbuff = new CarmaObj<T>(this->current_context, dims_data1);
  dims_data1[1] = this->nslope();
  this->d_compbuff2 = new CarmaObj<T>(this->current_context, dims_data1);
  long dims_data2[3] = {2, this->nactu(), nmodes};
  this->d_M2V = new CarmaObj<T>(this->current_context, dims_data2, M2V);
  dims_data2[1] = this->nslope();
  dims_data2[2] = nrec;
  this->d_slpol = new CarmaObj<T>(this->current_context, dims_data2);
  dims_data2[1] = nmodes;
  dims_data2[2] = this->nslope();
  this->d_S2M = new CarmaObj<T>(this->current_context, dims_data2);
  dims_data2[1] = this->nslope();
  dims_data2[2] = nmodes;
  CarmaObj<T> *d_tmp = new CarmaObj<T>(this->current_context, dims_data2);
  dims_data2[1] = nmodes;
  CarmaObj<T> *d_tmp2 = new CarmaObj<T>(this->current_context, dims_data2);

  std::cout << "Computing S2M matrix..." << std::endl;
  // 1. tmp = D*M2V
  carma_gemm(this->cublas_handle(), 'n', 'n', this->nslope(), nmodes,
             this->nactu(), 1.0f, this->d_imat->get_data(), this->nslope(),
             d_M2V->get_data(), this->nactu(), 0.0f, d_tmp->get_data(),
             this->nslope());
  // 2. tmp2 = (D*M2V)t * (D*M2V)
  carma_gemm(this->cublas_handle(), 't', 'n', nmodes, nmodes, this->nslope(),
             1.0f, d_tmp->get_data(), this->nslope(), d_tmp->get_data(),
             this->nslope(), 0.0f, d_tmp2->get_data(), nmodes);

  // 3. tmp2 = (tmp2)⁻¹
  carma_potr_inv(d_tmp2);
  // carma_potr_inv(d_tmp2);
  // 4. S2M = (D*M2V)⁻¹
  carma_gemm(this->cublas_handle(), 'n', 't', nmodes, this->nslope(), nmodes,
             1.0f, d_tmp2->get_data(), nmodes, d_tmp->get_data(), this->nslope(),
             0.0f, d_S2M->get_data(), nmodes);

  delete d_tmp;
  delete d_tmp2;

  std::cout << "Computing transfer functions..." << std::endl;
  compute_Hcor();

  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int sutra_controller_ls<T, Tout>::modalControlOptimization() {
  this->current_context->set_active_device(this->device, 1);
  if (this->cpt_rec >= this->delay) {
    // POLC to retrieve open-loop measurements for further refreshing modal
    // gains
    this->d_compbuff->copy(this->d_com_clipped, 1, 1);
    carma_gemv<T>(this->cublas_handle(), 'n', this->nslope(), this->nactu(),
                  1.0f, *d_imat, this->nslope(), *d_compbuff, 1, 0.0f,
                  *d_compbuff2, 1);
    carma_geam<T>(this->cublas_handle(), 'n', 'n', this->nslope(), 1, 1.0f,
                  *this->d_centroids, this->nslope(), -1.0f, *d_compbuff2,
                  this->nslope(),
                  this->d_slpol->get_data_at((this->cpt_rec - (int)this->delay) *
                                           this->nslope()),
                  this->nslope());
  }
  this->cpt_rec++;

  long dims_data[2] = {1, this->nrec / 2 + 1};
  CarmaObj<cuFloatComplex> d_FFT(this->current_context, dims_data);
  dims_data[1] = this->nrec / 2;
  CarmaObj<T> d_fftmodes(this->current_context, dims_data);
  dims_data[1] = this->ngain;
  CarmaObj<T> d_phaseError(this->current_context, dims_data);
  long dims_data2[3] = {2, this->nrec, this->nmodes};
  CarmaObj<T> d_modes(this->current_context, dims_data2);
  int imin;
  T mgain[this->nmodes];

  // 1. modes = S2M * slopes_open_loop and transpose for fft
  carma_gemm(this->cublas_handle(), 't', 't', this->nrec, this->nmodes,
             this->nslope(), 1.0f, this->d_slpol->get_data(), this->nslope(),
             this->d_S2M->get_data(), this->nmodes, 0.0f, d_modes.get_data(),
             this->nrec);
  this->d_slpol->scale(0.0f, 1);

  // 2. Init and compute FFT modes
  dims_data[1] = this->nrec;
  carma_initfft<T, cuFloatComplex>(dims_data, d_modes.get_plan(), CUFFT_R2C);
  for (int i = 0; i < this->nmodes; i++) {
    CarmaFFT<T, cuFloatComplex>(d_modes.get_data_at(i * this->nrec),
                                 d_FFT.get_data(), 1, *d_modes.get_plan());
    absnormfft(d_FFT.get_data(), d_fftmodes.get_data(), this->nrec / 2,
               2.0f / (T)this->nrec,
               this->current_context->get_device(this->device));
    carma_gemv(this->cublas_handle(), 'n', this->ngain, this->nrec / 2, 1.0f,
               this->d_Hcor->get_data(), this->ngain, d_fftmodes.get_data(), 1,
               0.0f, d_phaseError.get_data(), 1);

    // Find and store optimum gain for mode i
    imin = carma_where_amin(this->cublas_handle(), this->ngain,
                            d_phaseError.get_data(), 1) -
           1;
    mgain[i] =
        this->gmin + imin * (this->gmax - this->gmin) / (this->ngain - 1);
  }

  this->d_gain->host2device(mgain);
  // Compute CMAT
  build_cmat_modopti();

  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int sutra_controller_ls<T, Tout>::loadopen_loopSlp(T *ol_slopes) {
  this->current_context->set_active_device(this->device, 1);
  this->d_slpol->host2device(ol_slopes);

  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int sutra_controller_ls<T, Tout>::compute_Hcor() {
  this->current_context->set_active_device(this->device, 1);
  long dims_data[3] = {2, this->ngain, this->nrec / 2};
  this->d_Hcor = new CarmaObj<T>(this->current_context, dims_data);

  compute_Hcor_gpu(this->d_Hcor->get_data(), this->ngain, this->nrec / 2,
                   this->Fs, this->gmin, this->gmax, this->delay,
                   this->current_context->get_device(this->device));

  return EXIT_SUCCESS;
}

template class sutra_controller_ls<float, float>;
template class sutra_controller_ls<float, uint16_t>;
