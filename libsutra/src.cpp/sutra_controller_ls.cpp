#include <sutra_controller_ls.h>
#include <string>

sutra_controller_ls::sutra_controller_ls(carma_context *context, long nvalid,
                                         long nslope, long nactu, float delay,
                                         sutra_dms *dms, int *idx_dms, int ndm)
    : sutra_controller(context, nvalid, nslope, nactu, delay, dms, idx_dms,
                       ndm) {
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
  this->d_imat = new carma_obj<float>(context, dims_data2);

  dims_data2[1] = nactu;
  dims_data2[2] = nslope;
  this->d_cmat = new carma_obj<float>(context, dims_data2);

  dims_data2[1] = dims_data2[2] = nactu;
  d_U = new carma_obj<float>(current_context, dims_data2);

  dims_data1[1] = nslope < nactu ? nslope : nactu;
  this->d_eigenvals = new carma_obj<float>(context, dims_data1);
  this->h_eigenvals = new carma_host_obj<float>(dims_data1, MA_PAGELOCK);

  if (delay > 0) {
    dims_data2[1] = nslope;
    dims_data2[2] = (int)delay + 1;
    this->d_cenbuff = new carma_obj<float>(context, dims_data2);
  } else
    this->d_cenbuff = 0L;
  this->d_M2V = 0L;
  this->d_S2M = 0L;
  this->d_slpol = 0L;
  this->d_Hcor = 0L;
  this->d_com1 = 0L;
  this->d_com2 = 0L;
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
  this->d_err = new carma_obj<float>(context, dims_data1);
  this->d_gain = new carma_obj<float>(context, dims_data1);
}

sutra_controller_ls::~sutra_controller_ls() {
  current_context->set_activeDevice(device, 1);

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
    delete this->d_com1;
    delete this->d_com2;
    delete this->d_compbuff;
    delete this->d_compbuff2;
  }
}

string sutra_controller_ls::get_type() { return "ls"; }

int sutra_controller_ls::svdec_imat() {
  // doing U = Dt.D where D is i_mat
  float one = 1., zero = 0.;
  int nCols = this->d_imat->getDims(2);

  current_context->set_activeDevice(device, 1);
  if (d_U->getDims(1) !=
      nCols) {  // d_imat shape was changed during modal basis. Adapt yourself.
    delete d_U;
    this->d_U = new carma_obj<float>(current_context,
                                     std::vector<long>{2, nCols, nCols}.data());
    delete d_eigenvals;
    delete h_eigenvals;
    this->d_eigenvals = new carma_obj<float>(
        current_context, std::vector<long>{1, nCols}.data());
    this->h_eigenvals = new carma_host_obj<float>(
        std::vector<long>{1, nCols}.data(), MA_PAGELOCK);
  }

  if (carma_syrk<float>(cublas_handle(), CUBLAS_FILL_MODE_LOWER, 't', nCols,
                        nslope(), one, *d_imat, nslope(), zero, *d_U, nCols)) {
    return EXIT_FAILURE;
  }

  if (!magma_disabled()) {
    // we can skip this step syevd use only the lower part
    fill_sym_matrix('U', d_U->getData(), nCols, nCols * nCols,
                    current_context->get_device(device));

    // doing evd of U inplace
    if (carma_syevd<float, 1>('V', d_U, h_eigenvals) == EXIT_FAILURE) {
      // if (syevd_f('V', d_U, h_eigenvals) == EXIT_FAILURE) {
      // Case where MAGMA is not feeling good :-/
      return EXIT_FAILURE;
    }
    d_eigenvals->host2device(h_eigenvals->getData());
  } else {  // CULA case
    // We fill the upper matrix part of the matrix
    fill_sym_matrix<float>('L', *d_U, nCols, nCols * nCols,
                           current_context->get_device(device));

    carma_obj<float> d_tmp(d_U);
    carma_obj<float> d_tmp2(d_U);

    if (carma_cula_svd<float>(&d_tmp, d_eigenvals, d_U, &d_tmp2) ==
        EXIT_FAILURE) {
      return EXIT_FAILURE;
    }
    d_eigenvals->device2host(*h_eigenvals);
  }
  return EXIT_SUCCESS;
}

int sutra_controller_ls::set_gain(float gain) {
  this->gain = gain;
  return EXIT_SUCCESS;
}

int sutra_controller_ls::set_cmat(float *cmat) {
  current_context->set_activeDevice(device, 1);
  this->d_cmat->host2device(cmat);
  return EXIT_SUCCESS;
}

int sutra_controller_ls::set_imat(float *imat) {
  current_context->set_activeDevice(device, 1);
  this->d_imat->host2device(imat);
  return EXIT_SUCCESS;
}

int sutra_controller_ls::set_mgain(float *mgain) {
  current_context->set_activeDevice(device, 1);
  this->d_gain->host2device(mgain);
  return EXIT_SUCCESS;
}

int sutra_controller_ls::set_delay(float delay) {
  this->delay = delay;
  return EXIT_SUCCESS;
}

int sutra_controller_ls::build_cmat(int nfilt, bool filt_tt) {
  current_context->set_activeDevice(device, 1);

  long dims_data1[2] = {1, 0};
  long dims_data2[3] = {2, 0, 0};

  dims_data2[1] = dims_data2[2] = nactu();
  carma_obj<float> d_tmp(current_context, dims_data2),
      d_tmp2(current_context, dims_data2);

  dims_data1[1] = nactu();
  carma_obj<float> d_eigenvals_inv(current_context, dims_data1);
  carma_host_obj<float> h_eigenvals_inv(dims_data1, MA_PAGELOCK);

  float one = 1., zero = 0.;

  int nb_elem = this->h_eigenvals->getNbElem();
  memset(h_eigenvals_inv.getData(), 0, sizeof(float) * nb_elem);

  // filtering modes
  /*
  if (filt_tt)
    nb_elem -= 2;
  */
  if (!magma_disabled()) {
    for (int cc = nfilt; cc < nb_elem; cc++) {
      float eigenval = (*this->h_eigenvals)[cc];
      h_eigenvals_inv[cc] =  // 1.0f / eigenval;
          (fabs(eigenval) > 1.e-9) ? 1.0f / eigenval : 0.f;
    }
  } else {
    for (int cc = 0; cc < nb_elem - nfilt; cc++) {
      float eigenval = (*this->h_eigenvals)[cc];
      h_eigenvals_inv[cc] = (fabs(eigenval) > 1.e-9) ? 1.0f / eigenval : 0.f;
    }
  }
  d_eigenvals_inv.host2device(h_eigenvals_inv.getData());

  carma_dgmm(cublas_handle(), CUBLAS_SIDE_RIGHT, nactu(), nactu(),
             d_U->getData(), nactu(), d_eigenvals_inv.getData(), one,
             d_tmp.getData(), nactu());
  carma_gemm(cublas_handle(), 'n', 't', nactu(), nactu(), nactu(), one,
             d_tmp.getData(), nactu(), d_U->getData(), nactu(), zero,
             d_tmp2.getData(), nactu());
  carma_gemm(cublas_handle(), 'n', 't', nactu(), nslope(), nactu(), one,
             d_tmp2.getData(), nactu(), d_imat->getData(), nslope(), zero,
             d_cmat->getData(), nactu());

  return EXIT_SUCCESS;
}

int sutra_controller_ls::build_cmat(int nfilt) {
  current_context->set_activeDevice(device, 1);
  return this->build_cmat(nfilt, false);
}

int sutra_controller_ls::frame_delay() {
  // here we place the content of d_centroids into cenbuf and get
  // the actual centroid frame for error computation depending on delay value

  current_context->set_activeDevice(device, 1);
  if (delay > 0) {
    for (int cc = 0; cc < delay; cc++)
      shift_buf(this->d_cenbuff->getDataAt(cc * this->nslope()), 1,
                this->nslope(), this->current_context->get_device(device));

    carmaSafeCall(
        cudaMemcpy(this->d_cenbuff->getDataAt((int)delay * this->nslope()),
                   this->d_centroids->getData(), sizeof(float) * this->nslope(),
                   cudaMemcpyDeviceToDevice));

    carmaSafeCall(
        cudaMemcpy(this->d_centroids->getData(), this->d_cenbuff->getData(),
                   sizeof(float) * this->nslope(), cudaMemcpyDeviceToDevice));
  }

  return EXIT_SUCCESS;
}

int sutra_controller_ls::comp_com() {
  current_context->set_activeDevice(device, 1);

  // this->frame_delay();
  int nstreams = streams->get_nbStreams();

  // Modal Control Optimization
  if (this->is_modopti) {
    // Refresh when enough slopes have been recorded
    if (this->cpt_rec >= this->nrec + this->delay) {
      std::cout << "Refreshing modal gains..." << std::endl;
      modalControlOptimization();
      this->cpt_rec = 0;
    }
    if (cpt_rec >= this->delay) {
      // POLC to retrieve open-loop measurements for further refreshing modal
      // gains
      this->d_com2->copy(this->d_com1, 1, 1);
      this->d_com1->copy(this->d_com, 1, 1);
      // POLC equations
      carma_geam<float>(cublas_handle(), 'n', 'n', nactu(), 1,
                        (float)(delay - 1), this->d_com2->getData(), nactu(),
                        1.0f - (delay - 1), this->d_com1->getData(), nactu(),
                        this->d_compbuff->getData(), nactu());
      carma_gemv<float>(cublas_handle(), 'n', nslope(), nactu(), 1.0f, *d_imat,
                        nslope(), *d_compbuff, 1, 0.0f, *d_compbuff2, 1);
      carma_geam<float>(cublas_handle(), 'n', 'n', nslope(), 1, 1.0f,
                        *d_centroids, nslope(), -1.0f, *d_compbuff2, nslope(),
                        this->d_slpol->getDataAt(
                            (this->cpt_rec - (int)this->delay) * nslope()),
                        nslope());
    }
    this->cpt_rec++;
  }

  // INTEGRATOR
  if (nstreams > 1) {
    float alpha = -1.0f;
    float beta = 0.0f;

    for (int i = 0; i < nstreams; i++) {
      int istart1 =
          i * this->d_cmat->getDims(2) * this->d_cmat->getDims(1) / nstreams;
      int istart2 = i * this->d_cmat->getDims(1) / nstreams;

      // cout << istart1 << " " << istart2 << endl;

      cublasSetStream(cublas_handle(), this->streams->get_stream(i));

      cublasOperation_t trans = carma_char2cublasOperation('n');

      carma_checkCublasStatus(cublasSgemv(
          cublas_handle(), trans, this->d_cmat->getDims(1) / nstreams,
          this->d_cmat->getDims(2), &alpha,
          &((this->d_cmat->getData())[istart1]),
          this->d_cmat->getDims(1) / nstreams, this->d_centroids->getData(), 1,
          &beta, &((this->d_err->getData())[istart2]), 1));
    }

    mult_int(this->d_com->getData(), this->d_err->getData(),
             this->d_gain->getData(), this->gain, this->nactu(),
             this->current_context->get_device(device), this->streams);

    this->streams->wait_all_streams();

  } else {
    //    float *cmat=(float*)malloc(this->d_cmat->getNbElem()*sizeof(float));
    //    d_cmat->device2host(cmat);
    //    DEBUG_TRACE("here %f %f %d", cmat[0], this->gain, this->open_loop);
    // compute error
    this->d_err->gemv('n', -1.0f, this->d_cmat, this->d_cmat->getDims(1),
                      this->d_centroids, 1, 0.0f, 1);

    // apply modal gain & loop gain
    if (this->is_modopti)
      mult_int(this->d_com->getData(), this->d_err->getData(), this->gain,
               this->nactu(), this->current_context->get_device(device));
    else
      mult_int(this->d_com->getData(), this->d_err->getData(),
               this->d_gain->getData(), this->gain, this->nactu(),
               this->current_context->get_device(device));
  }

  return EXIT_SUCCESS;
}

int sutra_controller_ls::build_cmat_modopti() {
  current_context->set_activeDevice(device, 1);
  long dims_data2[3] = {2, nactu(), this->nmodes};
  carma_obj<float> d_tmp(current_context, dims_data2);

  // Compute cmat as M2V*(modal gains)*S2M
  carma_dgmm(cublas_handle(), CUBLAS_SIDE_RIGHT, nactu(), this->nmodes,
             this->d_M2V->getData(), nactu(), this->d_gain->getData(), 1,
             d_tmp.getData(), nactu());
  carma_gemm(cublas_handle(), 'n', 'n', nactu(), nslope(), this->nmodes, 1.0f,
             d_tmp.getData(), nactu(), this->d_S2M->getData(), this->nmodes,
             0.0f, this->d_cmat->getData(), nactu());

  return EXIT_SUCCESS;
}

int sutra_controller_ls::init_modalOpti(int nmodes, int nrec, float *M2V,
                                        float gmin, float gmax, int ngain,
                                        float Fs) {
  current_context->set_activeDevice(device, 1);
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
  this->d_gain = new carma_obj<float>(current_context, dims_data1);
  dims_data1[1] = nactu();
  this->d_com1 = new carma_obj<float>(current_context, dims_data1);
  this->d_com2 = new carma_obj<float>(current_context, dims_data1);
  this->d_compbuff = new carma_obj<float>(this->current_context, dims_data1);
  dims_data1[1] = nslope();
  this->d_compbuff2 = new carma_obj<float>(this->current_context, dims_data1);
  long dims_data2[3] = {2, nactu(), nmodes};
  this->d_M2V = new carma_obj<float>(current_context, dims_data2, M2V);
  dims_data2[1] = nslope();
  dims_data2[2] = nrec;
  this->d_slpol = new carma_obj<float>(current_context, dims_data2);
  dims_data2[1] = nmodes;
  dims_data2[2] = nslope();
  this->d_S2M = new carma_obj<float>(current_context, dims_data2);
  dims_data2[1] = nslope();
  dims_data2[2] = nmodes;
  carma_obj<float> *d_tmp = new carma_obj<float>(current_context, dims_data2);
  dims_data2[1] = nmodes;
  carma_obj<float> *d_tmp2 = new carma_obj<float>(current_context, dims_data2);

  std::cout << "Computing S2M matrix..." << std::endl;
  // 1. tmp = D*M2V
  carma_gemm(cublas_handle(), 'n', 'n', nslope(), nmodes, nactu(), 1.0f,
             this->d_imat->getData(), nslope(), d_M2V->getData(), nactu(), 0.0f,
             d_tmp->getData(), nslope());
  // 2. tmp2 = (D*M2V)t * (D*M2V)
  carma_gemm(cublas_handle(), 't', 'n', nmodes, nmodes, nslope(), 1.0f,
             d_tmp->getData(), nslope(), d_tmp->getData(), nslope(), 0.0f,
             d_tmp2->getData(), nmodes);

  // 3. tmp2 = (tmp2)⁻¹
  carma_potri(d_tmp2);
  // 4. S2M = (D*M2V)⁻¹
  carma_gemm(cublas_handle(), 'n', 't', nmodes, nslope(), nmodes, 1.0f,
             d_tmp2->getData(), nmodes, d_tmp->getData(), nslope(), 0.0f,
             d_S2M->getData(), nmodes);

  delete d_tmp;
  delete d_tmp2;

  std::cout << "Computing transfer functions..." << std::endl;
  compute_Hcor();

  return EXIT_SUCCESS;
}

int sutra_controller_ls::modalControlOptimization() {
  current_context->set_activeDevice(device, 1);
  long dims_data[2] = {1, this->nrec / 2 + 1};
  carma_obj<cuFloatComplex> d_FFT(current_context, dims_data);
  dims_data[1] = this->nrec / 2;
  carma_obj<float> d_fftmodes(current_context, dims_data);
  dims_data[1] = this->ngain;
  carma_obj<float> d_phaseError(current_context, dims_data);
  long dims_data2[3] = {2, this->nrec, this->nmodes};
  carma_obj<float> d_modes(current_context, dims_data2);
  int imin;
  float mgain[this->nmodes];

  // 1. modes = S2M * slopes_open_loop and transpose for fft
  carma_gemm(cublas_handle(), 't', 't', this->nrec, this->nmodes, nslope(),
             1.0f, this->d_slpol->getData(), nslope(), this->d_S2M->getData(),
             this->nmodes, 0.0f, d_modes.getData(), this->nrec);
  this->d_slpol->scale(0.0f, 1);

  // 2. Init and compute FFT modes
  dims_data[1] = this->nrec;
  carma_initfft<float, cuFloatComplex>(dims_data, d_modes.getPlan(), CUFFT_R2C);
  for (int i = 0; i < this->nmodes; i++) {
    carma_fft<float, cuFloatComplex>(d_modes.getDataAt(i * this->nrec),
                                     d_FFT.getData(), 1, *d_modes.getPlan());
    absnormfft(d_FFT.getData(), d_fftmodes.getData(), this->nrec / 2,
               2.0f / (float)this->nrec,
               this->current_context->get_device(device));
    carma_gemv(cublas_handle(), 'n', this->ngain, this->nrec / 2, 1.0f,
               this->d_Hcor->getData(), this->ngain, d_fftmodes.getData(), 1,
               0.0f, d_phaseError.getData(), 1);

    // Find and store optimum gain for mode i
    imin = carma_where_amin(cublas_handle(), this->ngain,
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

int sutra_controller_ls::loadOpenLoopSlp(float *ol_slopes) {
  current_context->set_activeDevice(device, 1);
  this->d_slpol->host2device(ol_slopes);

  return EXIT_SUCCESS;
}

int sutra_controller_ls::compute_Hcor() {
  current_context->set_activeDevice(device, 1);
  long dims_data[3] = {2, this->ngain, this->nrec / 2};
  this->d_Hcor = new carma_obj<float>(current_context, dims_data);

  compute_Hcor_gpu(this->d_Hcor->getData(), this->ngain, this->nrec / 2,
                   this->Fs, this->gmin, this->gmax, this->delay,
                   this->current_context->get_device(device));

  return EXIT_SUCCESS;
}
