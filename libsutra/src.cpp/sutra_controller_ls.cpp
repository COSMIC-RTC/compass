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

//! \file      sutra_controller_ls.cpp
//! \ingroup   libsutra
//! \class     sutra_controller_ls
//! \brief     this class provides the controller_ls features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.3.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#include <carma_magma.h>
#include <sutra_controller_ls.h>
#include <string>

template <typename T, typename Tout>
sutra_controller_ls<T, Tout>::sutra_controller_ls(carma_context *context,
                                                  long nvalid, long nslope,
                                                  long nactu, float delay,
                                                  sutra_dms *dms, int *idx_dms,
                                                  int ndm, int *idx_centro, int ncentro)
    : sutra_controller<T, Tout>(context, nvalid, nslope, nactu, delay, dms,
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
  this->d_imat = new carma_obj<T>(context, dims_data2);

  dims_data2[1] = nactu;
  dims_data2[2] = nslope;
  this->d_cmat = new carma_obj<T>(context, dims_data2);

  dims_data2[1] = dims_data2[2] = nactu;
  d_U = new carma_obj<T>(context, dims_data2);

  dims_data1[1] = nslope < nactu ? nslope : nactu;
  this->d_eigenvals = new carma_obj<T>(context, dims_data1);
  this->h_eigenvals = new carma_host_obj<T>(dims_data1, MA_PAGELOCK);

  if (delay > 0) {
    dims_data2[1] = nslope;
    dims_data2[2] = (int)delay + 1;
    this->d_cenbuff = new carma_obj<T>(context, dims_data2);
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
  this->d_err = new carma_obj<T>(context, dims_data1);
  this->d_gain = new carma_obj<T>(context, dims_data1);
}

template <typename T, typename Tout>
sutra_controller_ls<T, Tout>::~sutra_controller_ls() {
  this->current_context->set_activeDevice(this->device, 1);

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
  int nCols = this->d_imat->getDims(2);

  this->current_context->set_activeDevice(this->device, 1);
  if (d_U->getDims(1) !=
      nCols) {  // d_imat shape was changed during modal basis. Adapt yourself.
    delete d_U;
    this->d_U = new carma_obj<T>(this->current_context,
                                 std::vector<long>{2, nCols, nCols}.data());
    delete d_eigenvals;
    delete h_eigenvals;
    this->d_eigenvals = new carma_obj<T>(this->current_context,
                                         std::vector<long>{1, nCols}.data());
    this->h_eigenvals =
        new carma_host_obj<T>(std::vector<long>{1, nCols}.data(), MA_PAGELOCK);
  }

  if (carma_syrk<T>(this->cublas_handle(), CUBLAS_FILL_MODE_LOWER, 't', nCols,
                    this->nslope(), one, *d_imat, this->nslope(), zero, *d_U,
                    nCols)) {
    return EXIT_FAILURE;
  }

  if (!carma_magma_disabled()) {
    // we can skip this step syevd use only the lower part
    fill_sym_matrix('U', d_U->getData(), nCols, nCols * nCols,
                    this->current_context->get_device(this->device));

    // doing evd of U inplace
    if (carma_magma_syevd<T>('V', d_U, h_eigenvals) == EXIT_FAILURE) {
      // if (syevd_f('V', d_U, h_eigenvals) == EXIT_FAILURE) {
      // Case where MAGMA is not feeling good :-/
      return EXIT_FAILURE;
    }
    d_eigenvals->host2device(h_eigenvals->getData());
  } else {  // CULA case
    // We fill the upper matrix part of the matrix
    fill_sym_matrix<T>('L', *d_U, nCols, nCols * nCols,
                       this->current_context->get_device(this->device));

    carma_obj<T> d_tmp(d_U);
    carma_obj<T> d_tmp2(d_U);

    if (carma_cula_svd<T>(&d_tmp, d_eigenvals, d_U, &d_tmp2) == EXIT_FAILURE) {
      return EXIT_FAILURE;
    }
    d_eigenvals->device2host(h_eigenvals->getData());
  }
  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int sutra_controller_ls<T, Tout>::set_cmat(T *cmat) {
  this->current_context->set_activeDevice(this->device, 1);
  this->d_cmat->host2device(cmat);
  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int sutra_controller_ls<T, Tout>::set_imat(T *imat) {
  this->current_context->set_activeDevice(this->device, 1);
  this->d_imat->host2device(imat);
  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int sutra_controller_ls<T, Tout>::set_mgain(T *mgain) {
  this->current_context->set_activeDevice(this->device, 1);
  this->d_gain->host2device(mgain);
  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int sutra_controller_ls<T, Tout>::build_cmat(int nfilt, bool filt_tt) {
  this->current_context->set_activeDevice(this->device, 1);

  long dims_data1[2] = {1, 0};
  long dims_data2[3] = {2, 0, 0};

  dims_data2[1] = dims_data2[2] = this->nactu();
  carma_obj<T> d_tmp(this->current_context, dims_data2),
      d_tmp2(this->current_context, dims_data2);

  dims_data1[1] = this->nactu();
  carma_obj<T> d_eigenvals_inv(this->current_context, dims_data1);
  carma_host_obj<T> h_eigenvals_inv(dims_data1, MA_PAGELOCK);

  T one = 1., zero = 0.;

  int nb_elem = this->h_eigenvals->getNbElem();
  memset(h_eigenvals_inv.getData(), 0, sizeof(T) * nb_elem);

  // filtering modes
  /*
  if (filt_tt)
    nb_elem -= 2;
  */
  if (!carma_magma_disabled()) {
    for (int cc = nfilt; cc < nb_elem; cc++) {
      T eigenval = (*this->h_eigenvals)[cc];
      h_eigenvals_inv[cc] =  // 1.0f / eigenval;
          (fabs(eigenval) > 1.e-9) ? 1.0f / eigenval : 0.f;
    }
  } else {
    for (int cc = 0; cc < nb_elem - nfilt; cc++) {
      T eigenval = (*this->h_eigenvals)[cc];
      h_eigenvals_inv[cc] = (fabs(eigenval) > 1.e-9) ? 1.0f / eigenval : 0.f;
    }
  }
  d_eigenvals_inv.host2device(h_eigenvals_inv.getData());

  carma_dgmm(this->cublas_handle(), CUBLAS_SIDE_RIGHT, this->nactu(),
             this->nactu(), d_U->getData(), this->nactu(),
             d_eigenvals_inv.getData(), one, d_tmp.getData(), this->nactu());
  carma_gemm(this->cublas_handle(), 'n', 't', this->nactu(), this->nactu(),
             this->nactu(), one, d_tmp.getData(), this->nactu(), d_U->getData(),
             this->nactu(), zero, d_tmp2.getData(), this->nactu());
  carma_gemm(this->cublas_handle(), 'n', 't', this->nactu(), this->nslope(),
             this->nactu(), one, d_tmp2.getData(), this->nactu(),
             d_imat->getData(), this->nslope(), zero, d_cmat->getData(),
             this->nactu());

  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int sutra_controller_ls<T, Tout>::build_cmat(int nfilt) {
  this->current_context->set_activeDevice(this->device, 1);
  return this->build_cmat(nfilt, false);
}

template <typename T, typename Tout>
int sutra_controller_ls<T, Tout>::frame_delay() {
  // here we place the content of d_centroids into cenbuf and get
  // the actual centroid frame for error computation depending on delay value

  this->current_context->set_activeDevice(this->device, 1);
  if (this->delay > 0) {
    for (int cc = 0; cc < this->delay; cc++)
      shift_buf(this->d_cenbuff->getDataAt(cc * this->nslope()), 1,
                this->nslope(),
                this->current_context->get_device(this->device));

    carmaSafeCall(cudaMemcpy(
        this->d_cenbuff->getDataAt((int)this->delay * this->nslope()),
        this->d_centroids->getData(), sizeof(T) * this->nslope(),
        cudaMemcpyDeviceToDevice));

    carmaSafeCall(
        cudaMemcpy(this->d_centroids->getData(), this->d_cenbuff->getData(),
                   sizeof(T) * this->nslope(), cudaMemcpyDeviceToDevice));
  }

  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int sutra_controller_ls<T, Tout>::comp_com() {
  this->current_context->set_activeDevice(this->device, 1);

  // this->frame_delay();
  int nstreams = this->streams->get_nbStreams();

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
  if (nstreams > 1) {
    T alpha = -1.0f;
    T beta = 0.0f;

    for (int i = 0; i < nstreams; i++) {
      int istart1 =
          i * this->d_cmat->getDims(2) * this->d_cmat->getDims(1) / nstreams;
      int istart2 = i * this->d_cmat->getDims(1) / nstreams;

      // cout << istart1 << " " << istart2 << endl;

      cublasSetStream(this->cublas_handle(), this->streams->get_stream(i));

      cublasOperation_t trans = carma_char2cublasOperation('n');

      carma_checkCublasStatus(cublasSgemv(
          this->cublas_handle(), trans, this->d_cmat->getDims(1) / nstreams,
          this->d_cmat->getDims(2), &alpha,
          &((this->d_cmat->getData())[istart1]),
          this->d_cmat->getDims(1) / nstreams, this->d_centroids->getData(), 1,
          &beta, &((this->d_err->getData())[istart2]), 1));
    }

    mult_int(this->d_com->getData(), this->d_err->getData(),
             this->d_gain->getData(), this->gain, this->nactu(),
             this->current_context->get_device(this->device), this->streams);

    this->streams->wait_all_streams();

  } else {
    //    T *cmat=(T*)malloc(this->d_cmat->getNbElem()*sizeof(T));
    //    d_cmat->device2host(cmat);
    //    DEBUG_TRACE("here %f %f %d", cmat[0], this->gain, this->open_loop);
    // compute error
    this->d_err->gemv('n', -1.0f, this->d_cmat, this->d_cmat->getDims(1),
                      this->d_centroids, 1, 0.0f, 1);

    // apply modal gain & loop gain
    if (this->is_modopti)
      mult_int(this->d_com->getData(), this->d_err->getData(), this->gain,
               this->nactu(), this->current_context->get_device(this->device));
    else
      mult_int(this->d_com->getData(), this->d_err->getData(),
               this->d_gain->getData(), this->gain, this->nactu(),
               this->current_context->get_device(this->device));
  }

  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int sutra_controller_ls<T, Tout>::build_cmat_modopti() {
  this->current_context->set_activeDevice(this->device, 1);
  long dims_data2[3] = {2, this->nactu(), this->nmodes};
  carma_obj<T> d_tmp(this->current_context, dims_data2);

  // Compute cmat as M2V*(modal gains)*S2M
  carma_dgmm(this->cublas_handle(), CUBLAS_SIDE_RIGHT, this->nactu(),
             this->nmodes, this->d_M2V->getData(), this->nactu(),
             this->d_gain->getData(), 1, d_tmp.getData(), this->nactu());
  carma_gemm(this->cublas_handle(), 'n', 'n', this->nactu(), this->nslope(),
             this->nmodes, 1.0f, d_tmp.getData(), this->nactu(),
             this->d_S2M->getData(), this->nmodes, 0.0f,
             this->d_cmat->getData(), this->nactu());

  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int sutra_controller_ls<T, Tout>::init_modalOpti(int nmodes, int nrec, T *M2V,
                                                 T gmin, T gmax, int ngain,
                                                 T Fs) {
  this->current_context->set_activeDevice(this->device, 1);
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
  this->d_gain = new carma_obj<T>(this->current_context, dims_data1);
  dims_data1[1] = this->nactu();
  this->d_compbuff = new carma_obj<T>(this->current_context, dims_data1);
  dims_data1[1] = this->nslope();
  this->d_compbuff2 = new carma_obj<T>(this->current_context, dims_data1);
  long dims_data2[3] = {2, this->nactu(), nmodes};
  this->d_M2V = new carma_obj<T>(this->current_context, dims_data2, M2V);
  dims_data2[1] = this->nslope();
  dims_data2[2] = nrec;
  this->d_slpol = new carma_obj<T>(this->current_context, dims_data2);
  dims_data2[1] = nmodes;
  dims_data2[2] = this->nslope();
  this->d_S2M = new carma_obj<T>(this->current_context, dims_data2);
  dims_data2[1] = this->nslope();
  dims_data2[2] = nmodes;
  carma_obj<T> *d_tmp = new carma_obj<T>(this->current_context, dims_data2);
  dims_data2[1] = nmodes;
  carma_obj<T> *d_tmp2 = new carma_obj<T>(this->current_context, dims_data2);

  std::cout << "Computing S2M matrix..." << std::endl;
  // 1. tmp = D*M2V
  carma_gemm(this->cublas_handle(), 'n', 'n', this->nslope(), nmodes,
             this->nactu(), 1.0f, this->d_imat->getData(), this->nslope(),
             d_M2V->getData(), this->nactu(), 0.0f, d_tmp->getData(),
             this->nslope());
  // 2. tmp2 = (D*M2V)t * (D*M2V)
  carma_gemm(this->cublas_handle(), 't', 'n', nmodes, nmodes, this->nslope(),
             1.0f, d_tmp->getData(), this->nslope(), d_tmp->getData(),
             this->nslope(), 0.0f, d_tmp2->getData(), nmodes);

  // 3. tmp2 = (tmp2)⁻¹
  carma_magma_potri(d_tmp2);
  // 4. S2M = (D*M2V)⁻¹
  carma_gemm(this->cublas_handle(), 'n', 't', nmodes, this->nslope(), nmodes,
             1.0f, d_tmp2->getData(), nmodes, d_tmp->getData(), this->nslope(),
             0.0f, d_S2M->getData(), nmodes);

  delete d_tmp;
  delete d_tmp2;

  std::cout << "Computing transfer functions..." << std::endl;
  compute_Hcor();

  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int sutra_controller_ls<T, Tout>::modalControlOptimization() {
  this->current_context->set_activeDevice(this->device, 1);
  if (this->cpt_rec >= this->delay) {
    // POLC to retrieve open-loop measurements for further refreshing modal
    // gains
    carma_geam<T>(this->cublas_handle(), 'n', 'n', this->nactu(), 1,
                  (T)(this->delay - 1), this->d_com1->getData(), this->nactu(),
                  1.0f - (this->delay - 1), this->d_com->getData(),
                  this->nactu(), this->d_compbuff->getData(), this->nactu());
    carma_gemv<T>(this->cublas_handle(), 'n', this->nslope(), this->nactu(),
                  1.0f, *d_imat, this->nslope(), *d_compbuff, 1, 0.0f,
                  *d_compbuff2, 1);
    carma_geam<T>(this->cublas_handle(), 'n', 'n', this->nslope(), 1, 1.0f,
                  *this->d_centroids, this->nslope(), -1.0f, *d_compbuff2,
                  this->nslope(),
                  this->d_slpol->getDataAt((this->cpt_rec - (int)this->delay) *
                                           this->nslope()),
                  this->nslope());
  }
  this->cpt_rec++;

  long dims_data[2] = {1, this->nrec / 2 + 1};
  carma_obj<cuFloatComplex> d_FFT(this->current_context, dims_data);
  dims_data[1] = this->nrec / 2;
  carma_obj<T> d_fftmodes(this->current_context, dims_data);
  dims_data[1] = this->ngain;
  carma_obj<T> d_phaseError(this->current_context, dims_data);
  long dims_data2[3] = {2, this->nrec, this->nmodes};
  carma_obj<T> d_modes(this->current_context, dims_data2);
  int imin;
  T mgain[this->nmodes];

  // 1. modes = S2M * slopes_open_loop and transpose for fft
  carma_gemm(this->cublas_handle(), 't', 't', this->nrec, this->nmodes,
             this->nslope(), 1.0f, this->d_slpol->getData(), this->nslope(),
             this->d_S2M->getData(), this->nmodes, 0.0f, d_modes.getData(),
             this->nrec);
  this->d_slpol->scale(0.0f, 1);

  // 2. Init and compute FFT modes
  dims_data[1] = this->nrec;
  carma_initfft<T, cuFloatComplex>(dims_data, d_modes.getPlan(), CUFFT_R2C);
  for (int i = 0; i < this->nmodes; i++) {
    carma_fft<T, cuFloatComplex>(d_modes.getDataAt(i * this->nrec),
                                 d_FFT.getData(), 1, *d_modes.getPlan());
    absnormfft(d_FFT.getData(), d_fftmodes.getData(), this->nrec / 2,
               2.0f / (T)this->nrec,
               this->current_context->get_device(this->device));
    carma_gemv(this->cublas_handle(), 'n', this->ngain, this->nrec / 2, 1.0f,
               this->d_Hcor->getData(), this->ngain, d_fftmodes.getData(), 1,
               0.0f, d_phaseError.getData(), 1);

    // Find and store optimum gain for mode i
    imin = carma_where_amin(this->cublas_handle(), this->ngain,
                            d_phaseError.getData(), 1) -
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
int sutra_controller_ls<T, Tout>::loadOpenLoopSlp(T *ol_slopes) {
  this->current_context->set_activeDevice(this->device, 1);
  this->d_slpol->host2device(ol_slopes);

  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int sutra_controller_ls<T, Tout>::compute_Hcor() {
  this->current_context->set_activeDevice(this->device, 1);
  long dims_data[3] = {2, this->ngain, this->nrec / 2};
  this->d_Hcor = new carma_obj<T>(this->current_context, dims_data);

  compute_Hcor_gpu(this->d_Hcor->getData(), this->ngain, this->nrec / 2,
                   this->Fs, this->gmin, this->gmax, this->delay,
                   this->current_context->get_device(this->device));

  return EXIT_SUCCESS;
}

template class sutra_controller_ls<float, float>;
template class sutra_controller_ls<float, uint16_t>;
