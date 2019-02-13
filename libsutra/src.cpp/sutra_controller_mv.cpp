#include <sutra_controller_mv.h>
#include <sutra_controller_utils.h>
#include <string>

template <typename Tcomp, typename Tout>
sutra_controller_mv<Tcomp, Tout>::sutra_controller_mv(
    carma_context *context, long nvalid_, long nslope_, long nactu_,
    float delay, sutra_dms *dms, int *idx_dms, int ndm)
    : sutra_controller<Tcomp, Tout>(context, nvalid_, nslope_, nactu_, delay,
                                    dms, idx_dms, ndm) {
  this->gain = 0.0f;

  //  this->nstreams = 1; //nvalid/10;
  //  while (this->nactu() % this->nstreams != 0)
  //    nstreams--;
  //  std::cerr << "controller uses " << nstreams << " streams" << std::endl;
  //  streams = new carma_streams(nstreams);
  long dims_data2[3];
  dims_data2[0] = 2;
  dims_data2[1] = this->nslope();
  dims_data2[2] = this->nactu();
  this->d_imat = new carma_obj<Tcomp>(this->current_context, dims_data2);
  dims_data2[1] = this->nactu();
  dims_data2[2] = this->nslope();
  this->d_cmat = new carma_obj<Tcomp>(this->current_context, dims_data2);
  // dims_data2[1] = dims_data2[2] = this->nactu();
  // d_U = new carma_obj<Tcomp>(this->current_context, dims_data2);
  this->d_cenbuff = 0L;
  if ((int)delay > 0) {
    dims_data2[1] = this->nslope();
    dims_data2[2] = (int)delay + 1;
    this->d_cenbuff = new carma_obj<Tcomp>(this->current_context, dims_data2);
  }
  dims_data2[1] = this->nslope();
  dims_data2[2] = this->nslope();
  this->d_Cmm = new carma_obj<Tcomp>(this->current_context, dims_data2);

  long dims_data1[2];
  dims_data1[0] = 1;

  // dims_data1[1] = this->nslope() < this->nactu() ? this->nslope() :
  // this->nactu();
  this->h_eigenvals = 0L;
  this->h_Cmmeigenvals = 0L;

  dims_data1[1] = this->nslope();
  this->d_noisemat = new carma_obj<Tcomp>(this->current_context, dims_data1);
  this->d_olmeas = new carma_obj<Tcomp>(this->current_context, dims_data1);
  // this->d_compbuff = new carma_obj<Tcomp>(this->current_context, dims_data1);
  this->d_compbuff2 = new carma_obj<Tcomp>(this->current_context, dims_data1);
  dims_data1[1] = this->nactu();
  this->d_compbuff = new carma_obj<Tcomp>(this->current_context, dims_data1);
  this->d_err = new carma_obj<Tcomp>(this->current_context, dims_data1);
  this->d_gain = new carma_obj<Tcomp>(this->current_context, dims_data1);

  // Florian features
  dims_data2[1] = this->nactu();  // 64564;
  dims_data2[2] = this->nactu();
  this->d_covmat = new carma_obj<Tcomp>(this->current_context, dims_data2);
  dims_data2[1] = this->nactu();
  dims_data2[2] = this->nactu();
  this->d_KLbasis = new carma_obj<Tcomp>(this->current_context, dims_data2);
  this->d_Cphim = 0L;

  cublas_handle = this->current_context->get_cublasHandle();
  // carma_checkCublasStatus(cublasCreate(&(this->cublas_handle)));
}

template <typename Tcomp, typename Tout>
sutra_controller_mv<Tcomp, Tout>::~sutra_controller_mv() {
  this->current_context->set_activeDevice(this->device, 1);
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
int sutra_controller_mv<Tcomp, Tout>::set_gain(Tcomp gain) {
  this->gain = gain;
  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int sutra_controller_mv<Tcomp, Tout>::set_mgain(Tcomp *mgain) {
  this->current_context->set_activeDevice(this->device, 1);
  this->d_gain->host2device(mgain);
  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int sutra_controller_mv<Tcomp, Tout>::set_cmat(Tcomp *cmat) {
  this->current_context->set_activeDevice(this->device, 1);
  this->d_cmat->host2device(cmat);
  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int sutra_controller_mv<Tcomp, Tout>::set_imat(Tcomp *imat) {
  this->current_context->set_activeDevice(this->device, 1);
  this->d_imat->host2device(imat);
  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int sutra_controller_mv<Tcomp, Tout>::load_noisemat(Tcomp *noise) {
  this->current_context->set_activeDevice(this->device, 1);
  this->d_noisemat->host2device(noise);
  return EXIT_SUCCESS;
}
// Florian features
template <typename Tcomp, typename Tout>
int sutra_controller_mv<Tcomp, Tout>::compute_Cmm(
    sutra_atmos *atmos, sutra_sensors *sensors, double *L0, double *cn2,
    double *alphaX, double *alphaY, double diamTel, double cobs) {
  this->current_context->set_activeDevice(this->device, 1);

  struct gtomo_struct g_tomo;
  init_tomo_gpu_gb(&g_tomo, atmos, sensors, diamTel, cobs);
  update_tomo_sys_gpu_gb(&g_tomo, sensors, alphaX, alphaY);
  update_tomo_atm_gpu_gb(&g_tomo, sensors, atmos, L0, cn2, alphaX, alphaY);
  matcov_gpu_4(this->d_Cmm->getData(), this->nslope(), this->nslope(), 0, 0,
               this->nslope(), &g_tomo, atmos, sensors, alphaX, alphaY);
  free_tomo_gpu_gb(&g_tomo);

  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int sutra_controller_mv<Tcomp, Tout>::compute_Cphim(
    sutra_atmos *atmos, sutra_sensors *sensors, sutra_dms *dms, double *L0,
    double *cn2, double *alphaX, double *alphaY, double *X, double *Y,
    double *xactu, double *yactu, double diamTel, double *k2, long *NlayerDm,
    long *indLayerDm, double FoV, double *pitch, double *alt_dm) {
  this->current_context->set_activeDevice(this->device, 1);

  // Find number of actuators without TTcomp DM
  int Nactu = 0;
  vector<sutra_dm *>::iterator p;
  p = dms->d_dms.begin();
  while (p != dms->d_dms.end()) {
    sutra_dm *dm = *p;
    if (dm->type != "tt") {
      Nactu += dm->nactus;
    }
    p++;
  }

  long dims_data2[3] = {2, Nactu, this->nslope()};
  if (this->d_Cphim != 0L) delete this->d_Cphim;
  this->d_Cphim = new carma_obj<Tcomp>(this->current_context, dims_data2);

  struct cphim_struct cphim_struct;

  // Compute Cphim matrix
  init_cphim_struct(&cphim_struct, atmos, sensors, dms, diamTel);
  update_cphim_sys(&cphim_struct, sensors, alphaX, alphaY, xactu, yactu, X, Y,
                   NlayerDm, indLayerDm, alt_dm, pitch, k2, FoV);
  update_cphim_atm(&cphim_struct, sensors, atmos, L0, cn2, alphaX, alphaY);
  CPHIM(this->d_Cphim->getData(), Nactu, this->nslope(), 0, 0, Nactu,
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
  long dims_data2[3] = {2, this->d_Cphim->getDims()[1],
                        this->d_Cphim->getDims()[1]};
  carma_obj<Tcomp> *d_Nact =
      new carma_obj<Tcomp>(this->current_context, dims_data2);
  dims_data2[2] = this->nslope();
  carma_obj<Tcomp> *d_tmp =
      new carma_obj<Tcomp>(this->current_context, dims_data2);
  d_Nact->host2device(Nact);
  carma_potri(d_Nact);
  carma_gemm(cublas_handle, 'n', 'n', this->d_Cphim->getDims()[1],
             this->nslope(), this->d_Cphim->getDims()[1], 1.0f,
             d_Nact->getData(), this->d_Cphim->getDims()[1],
             this->d_Cphim->getData(), this->d_Cphim->getDims()[1], 0.0f,
             d_tmp->getData(), this->d_Cphim->getDims()[1]);
  this->d_Cphim->copy(d_tmp, 1, 1);

  delete d_Nact;
  delete d_tmp;

  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int sutra_controller_mv<Tcomp, Tout>::do_covmat(sutra_dm *ydm, char *method,
                                                int *indx_pup, long dim,
                                                Tcomp *xpos, Tcomp *ypos,
                                                long Nkl, Tcomp norm,
                                                Tcomp ampli) {
  this->current_context->set_activeDevice(this->device, 1);
  long dims_data[3];
  dims_data[0] = 2;
  dims_data[1] = Nkl;
  dims_data[2] = Nkl;
  carma_obj<Tcomp> *d_statcov =
      new carma_obj<Tcomp>(this->current_context, dims_data);
  long dims_data2[2];
  dims_data2[0] = 1;
  dims_data2[1] = this->nactu();
  carma_obj<Tcomp> *d_KLcov =
      new carma_obj<Tcomp>(this->current_context, dims_data2);

  dims_data2[1] = dim;
  carma_obj<int> *d_indx =
      new carma_obj<int>(this->current_context, dims_data2);

  // Compute the statistic matrix from actuators positions & Kolmogorov
  // statistic
  carma_obj<Tcomp> *d_xpos =
      new carma_obj<Tcomp>(this->current_context, dims_data);
  carma_obj<Tcomp> *d_ypos =
      new carma_obj<Tcomp>(this->current_context, dims_data);
  d_xpos->host2device(xpos);
  d_ypos->host2device(ypos);
  do_statmat(d_statcov->getData(), Nkl, d_xpos->getData(), d_ypos->getData(),
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
     carma_obj<Tcomp> *d_IF = new carma_obj<Tcomp>(this->current_context,
     dims_data); dims_data[1] = this->nactu(); dims_data[2] = this->nactu();
     carma_obj<Tcomp> *d_geocov = new carma_obj<Tcomp>(current_context,
     dims_data);

     // Get influence functions of the DM (to be CUsparsed)
     d_indx->host2device(indx_pup);
     ydm->get_IF(d_IF->getData(),d_indx->getData(),dim);

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
  carma_host_obj<Tcomp> *h_eigenvals =
      new carma_host_obj<Tcomp>(dims_data2, MA_PAGELOCK);

  carma_syevd<Tcomp>('N', d_statcov, h_eigenvals);

  if (ydm->type == "kl") {
    dims_data2[1] = this->nactu();
    carma_host_obj<Tcomp> h_KLcov(dims_data2, MA_PAGELOCK);
    if (strcmp(method, "inv") == 0) {
      for (int i = 0; i < this->nactu(); i++) {
        h_KLcov[i] = -1. / (h_eigenvals->getData())[i];
      }
      if (Nkl == this->nactu()) {
        h_KLcov[this->nactu() - 1] = 0.;
      }
      d_KLcov->host2device(h_KLcov.getData());
      add_md(this->d_covmat->getData(), this->d_covmat->getData(),
             d_KLcov->getData(), this->nactu(),
             this->current_context->get_device(this->device));
    }
    if (strcmp(method, "n") == 0) {
      for (int i = 0; i < this->nactu(); i++) {
        h_KLcov[i] = -(h_eigenvals->getData())[i];
      }
      d_KLcov->host2device(h_KLcov.getData());
      add_md(this->d_covmat->getData(), this->d_covmat->getData(),
             d_KLcov->getData(), this->nactu(),
             this->current_context->get_device(this->device));
    }
    delete h_KLcov;
  }
  if (ydm->type == "pzt") {
    if (strcmp(method, "inv") == 0) {
      // Inversion of the KL covariance matrix
      carma_host_obj<Tcomp> *h_eigenvals_inv =
          new carma_host_obj<Tcomp>(dims_data2, MA_PAGELOCK);
      for (int i = 0; i < this->nactu(); i++) {
        // Valeurs propres négatives.... A voir et inverser ordre si valeurs
        // propres positives
        h_eigenvals_inv->getData()[i] = -1. / h_eigenvals->getData()[i];
      }

      h_eigenvals_inv->getData()[this->nactu() - 1] = 0.;
      d_KLcov->host2device(h_eigenvals_inv->getData());

      // 2. Inversion of the KL basis
      dims_data2[0] = 1;
      dims_data2[1] = this->nactu();
      carma_obj<Tcomp> *d_eigen =
          new carma_obj<Tcomp>(this->current_context, dims_data2);
      dims_data[1] = this->nactu();
      dims_data[2] = this->nactu();
      carma_obj<Tcomp> *d_tmp =
          new carma_obj<Tcomp>(this->current_context, dims_data);
      carma_obj<Tcomp> *d_Ukl =
          new carma_obj<Tcomp>(this->current_context, dims_data);
      carma_obj<Tcomp> *d_Vkl =
          new carma_obj<Tcomp>(this->current_context, dims_data);
      carma_host_obj<Tcomp> *h_KL =
          new carma_host_obj<Tcomp>(dims_data, MA_PAGELOCK);
      carma_host_obj<Tcomp> *h_U =
          new carma_host_obj<Tcomp>(dims_data, MA_PAGELOCK);
      carma_host_obj<Tcomp> *h_Vt =
          new carma_host_obj<Tcomp>(dims_data, MA_PAGELOCK);

      h_KL->cpy_obj(this->d_KLbasis, cudaMemcpyDeviceToHost);

      carma_svd_cpu<Tcomp>(h_KL, h_eigenvals, h_U, h_Vt);

      d_Ukl->host2device(h_Vt->getData());
      d_Vkl->host2device(h_U->getData());

      for (int i = 0; i < this->nactu(); i++) {
        h_eigenvals_inv->getData()[i] = 1. / h_eigenvals->getData()[i];
      }

      d_eigen->host2device(h_eigenvals_inv->getData());

      carma_dgmm(cublas_handle, CUBLAS_SIDE_LEFT, this->nactu(), this->nactu(),
                 d_Vkl->getData(), this->nactu(), d_eigen->getData(), 1,
                 d_tmp->getData(), this->nactu());
      carma_gemm(cublas_handle, 't', 't', this->nactu(), this->nactu(),
                 this->nactu(), 1.0f, d_tmp->getData(), this->nactu(),
                 d_Ukl->getData(), this->nactu(), 0.0f,
                 this->d_KLbasis->getData(), this->nactu());

      // 3. Projection of KL covariance matrix in the DM basis

      carma_dgmm(cublas_handle, CUBLAS_SIDE_LEFT, this->nactu(), this->nactu(),
                 d_KLbasis->getData(), this->nactu(), d_KLcov->getData(), 1,
                 d_tmp->getData(), this->nactu());
      carma_gemm(cublas_handle, 't', 'n', this->nactu(), this->nactu(),
                 this->nactu(), 1.0f, d_tmp->getData(), this->nactu(),
                 this->d_KLbasis->getData(), this->nactu(), 0.0f,
                 this->d_covmat->getData(), this->nactu());

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
      carma_obj<Tcomp> *d_tmp =
          new carma_obj<Tcomp>(this->current_context, dims_data);
      for (int i = 0; i < this->nactu(); i++) {
        // Valeurs propres négatives.... A voir et inverser ordre si valeurs
        // propres positives
        h_eigenvals->getData()[i] = -h_eigenvals->getData()[i];
        std::cout << h_eigenvals->getData()[i] << std::endl;
      }
      h_eigenvals->getData()[this->nactu() - 1] = 0.;
      d_KLcov->host2device(h_eigenvals->getData());

      carma_dgmm(cublas_handle, CUBLAS_SIDE_RIGHT, this->nactu(), this->nactu(),
                 d_KLbasis->getData(), this->nactu(), d_KLcov->getData(), 1,
                 d_tmp->getData(), this->nactu());
      carma_gemm(cublas_handle, 'n', 't', this->nactu(), this->nactu(),
                 this->nactu(), 1.0f, d_tmp->getData(), this->nactu(),
                 this->d_KLbasis->getData(), this->nactu(), 0.0f,
                 this->d_covmat->getData(), this->nactu());

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
int sutra_controller_mv<Tcomp, Tout>::do_geomat(carma_obj<Tcomp> *d_geocov,
                                                carma_obj<Tcomp> *d_IF,
                                                long n_pts, Tcomp ampli) {
  this->current_context->set_activeDevice(this->device, 1);
  carma_gemm(cublas_handle, 't', 'n', this->nactu(), this->nactu(), n_pts, 1.0f,
             d_IF->getData(), n_pts, d_IF->getData(), n_pts, 0.0f,
             d_geocov->getData(), this->nactu());
  d_geocov->scale(ampli, 1);
  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int sutra_controller_mv<Tcomp, Tout>::piston_filt(carma_obj<Tcomp> *d_statcov) {
  this->current_context->set_activeDevice(this->device, 1);
  long Nmod = d_statcov->getDims()[1];
  long dims_data[3];
  dims_data[0] = 2;
  dims_data[1] = Nmod;
  dims_data[2] = Nmod;
  carma_obj<Tcomp> *d_F =
      new carma_obj<Tcomp>(this->current_context, dims_data);
  carma_obj<Tcomp> *d_tmp =
      new carma_obj<Tcomp>(this->current_context, dims_data);

  int N = d_statcov->getDims()[1] * d_statcov->getDims()[1];
  fill_filtmat(d_F->getData(), Nmod, N,
               this->current_context->get_device(this->device));

  carma_gemm(cublas_handle, 'n', 'n', Nmod, Nmod, Nmod, 1.0f, d_F->getData(),
             Nmod, d_statcov->getData(), Nmod, 0.0f, d_tmp->getData(), Nmod);
  carma_gemm(cublas_handle, 'n', 'n', Nmod, Nmod, Nmod, 1.0f, d_tmp->getData(),
             Nmod, d_F->getData(), Nmod, 0.0f, d_statcov->getData(), Nmod);

  delete d_tmp;
  delete d_F;

  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int sutra_controller_mv<Tcomp, Tout>::piston_filt_cphim(
    carma_obj<Tcomp> *d_cphim, Tcomp *F) {
  this->current_context->set_activeDevice(this->device, 1);

  long Nmod = d_cphim->getDims()[1];
  long dims_data[3];
  dims_data[0] = 2;
  dims_data[1] = Nmod;
  dims_data[2] = Nmod;
  carma_obj<Tcomp> *d_F =
      new carma_obj<Tcomp>(this->current_context, dims_data);
  dims_data[2] = this->nslope();
  carma_obj<Tcomp> *d_tmp =
      new carma_obj<Tcomp>(this->current_context, dims_data);

  d_F->host2device(F);

  carma_gemm(cublas_handle, 'n', 'n', Nmod, this->nslope(), Nmod, 1.0f,
             d_F->getData(), Nmod, d_cphim->getData(), Nmod, 0.0f,
             d_tmp->getData(), Nmod);
  d_cphim->copy(d_tmp, 1, 1);

  delete d_tmp;
  delete d_F;

  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int sutra_controller_mv<Tcomp, Tout>::invgen(carma_obj<Tcomp> *d_mat,
                                             Tcomp cond, int job) {
  this->current_context->set_activeDevice(this->device, 1);
  const long dims_data[3] = {2, d_mat->getDims()[1], d_mat->getDims()[2]};
  carma_obj<Tcomp> *d_U =
      new carma_obj<Tcomp>(this->current_context, dims_data);
  carma_obj<Tcomp> *d_tmp =
      new carma_obj<Tcomp>(this->current_context, dims_data);
  int i;

  const long dims_data2[2] = {1, d_mat->getDims()[1]};
  carma_obj<Tcomp> *d_eigenvals_inv =
      new carma_obj<Tcomp>(this->current_context, dims_data2);
  carma_host_obj<Tcomp> *h_eigenvals =
      new carma_host_obj<Tcomp>(dims_data2, MA_PAGELOCK);
  carma_host_obj<Tcomp> *h_eigenvals_inv =
      new carma_host_obj<Tcomp>(dims_data2, MA_PAGELOCK);

  d_U->copy(d_mat, 1, 1);
  carma_syevd<Tcomp>('V', d_U, h_eigenvals);
  // syevd_f('V',d_U,h_eigenvals);
  if (job == 1) {  // Conditionnement
    Tcomp maxe = h_eigenvals->getData()[d_mat->getDims()[1] - 1];

    for (i = 0; i < d_mat->getDims()[1]; i++) {
      if (h_eigenvals->getData()[i] < maxe / cond)
        h_eigenvals_inv->getData()[i] = 0.;
      else
        h_eigenvals_inv->getData()[i] = 1. / h_eigenvals->getData()[i];
    }
  }
  if (job == 0) {  // Filtre #cond modes
    for (i = 0; i < d_mat->getDims()[1]; i++) {
      if (i < cond)
        h_eigenvals_inv->getData()[i] = 0.;
      else
        h_eigenvals_inv->getData()[i] = 1. / h_eigenvals->getData()[i];
    }
  }

  d_eigenvals_inv->host2device(h_eigenvals_inv->getData());

  carma_dgmm(cublas_handle, CUBLAS_SIDE_RIGHT, d_mat->getDims()[1],
             d_mat->getDims()[2], d_U->getData(), d_mat->getDims()[1],
             d_eigenvals_inv->getData(), 1, d_tmp->getData(),
             d_mat->getDims()[1]);
  carma_gemm<Tcomp>(cublas_handle, 'n', 't', d_mat->getDims()[1],
                    d_mat->getDims()[1], d_mat->getDims()[2], 1.0f,
                    d_tmp->getData(), d_mat->getDims()[1], d_U->getData(),
                    d_mat->getDims()[1], 0.0f, d_mat->getData(),
                    d_mat->getDims()[1]);

  delete d_U;
  delete d_tmp;
  delete d_eigenvals_inv;
  delete h_eigenvals;
  delete h_eigenvals_inv;

  return EXIT_SUCCESS;
}
template <typename Tcomp, typename Tout>
int sutra_controller_mv<Tcomp, Tout>::invgen(carma_obj<Tcomp> *d_mat,
                                             carma_host_obj<Tcomp> *h_eigen,
                                             Tcomp cond) {
  this->current_context->set_activeDevice(this->device, 1);
  const long dims_data[3] = {2, d_mat->getDims()[1], d_mat->getDims()[2]};
  carma_obj<Tcomp> *d_U =
      new carma_obj<Tcomp>(this->current_context, dims_data);
  carma_obj<Tcomp> *d_tmp =
      new carma_obj<Tcomp>(this->current_context, dims_data);

  const long dims_data2[2] = {1, d_mat->getDims()[1]};
  carma_obj<Tcomp> *d_eigenvals_inv =
      new carma_obj<Tcomp>(this->current_context, dims_data2);
  carma_host_obj<Tcomp> *h_eigenvals_inv =
      new carma_host_obj<Tcomp>(dims_data2, MA_PAGELOCK);

  d_U->copy(d_mat, 1, 1);

  carma_syevd<Tcomp>('V', d_U, h_eigen);

  // syevd_f('V',d_U,h_eigen);
  // Conditionnement
  Tcomp maxe = h_eigen->getData()[d_mat->getDims()[1] - 3];
  int cpt = 0;
  for (int i = 0; i < d_mat->getDims()[1]; i++) {
    if (h_eigen->getData()[i] < maxe / cond) {
      h_eigenvals_inv->getData()[i] = 0.;
      cpt++;
    } else
      h_eigenvals_inv->getData()[i] = 1. / h_eigen->getData()[i];
  }

  d_eigenvals_inv->host2device(h_eigenvals_inv->getData());

  carma_dgmm(cublas_handle, CUBLAS_SIDE_RIGHT, d_mat->getDims()[1],
             d_mat->getDims()[2], d_U->getData(), d_mat->getDims()[1],
             d_eigenvals_inv->getData(), 1, d_tmp->getData(),
             d_mat->getDims()[1]);
  carma_gemm<Tcomp>(cublas_handle, 'n', 't', d_mat->getDims()[1],
                    d_mat->getDims()[1], d_mat->getDims()[2], 1.0f,
                    d_tmp->getData(), d_mat->getDims()[1], d_U->getData(),
                    d_mat->getDims()[1], 0.0f, d_mat->getData(),
                    d_mat->getDims()[1]);

  std::cout << "Inversion done with " << cpt << " modes filtered" << std::endl;

  delete d_U;
  delete d_tmp;
  delete d_eigenvals_inv;
  delete h_eigenvals_inv;

  return EXIT_SUCCESS;
}
template <typename Tcomp, typename Tout>
int sutra_controller_mv<Tcomp, Tout>::invgen_cpu(carma_obj<Tcomp> *d_mat,
                                                 carma_host_obj<Tcomp> *h_eigen,
                                                 Tcomp cond) {
  this->current_context->set_activeDevice(this->device, 1);
  const long dims_data[3] = {2, d_mat->getDims()[1], d_mat->getDims()[2]};
  carma_obj<Tcomp> *d_U =
      new carma_obj<Tcomp>(this->current_context, dims_data);
  carma_host_obj<Tcomp> *h_U =
      new carma_host_obj<Tcomp>(dims_data, MA_PAGELOCK);
  carma_host_obj<Tcomp> *h_V =
      new carma_host_obj<Tcomp>(dims_data, MA_PAGELOCK);
  carma_host_obj<Tcomp> *h_mat =
      new carma_host_obj<Tcomp>(dims_data, MA_PAGELOCK);
  carma_obj<Tcomp> *d_tmp =
      new carma_obj<Tcomp>(this->current_context, dims_data);

  const long dims_data2[2] = {1, d_mat->getDims()[1]};
  carma_obj<Tcomp> *d_eigenvals_inv =
      new carma_obj<Tcomp>(this->current_context, dims_data2);
  carma_host_obj<Tcomp> *h_eigenvals_inv =
      new carma_host_obj<Tcomp>(dims_data2, MA_PAGELOCK);

  d_mat->device2host(h_mat->getData());

  carma_svd_cpu<Tcomp>(h_mat, h_eigen, h_U, h_V);
  d_U->host2device(h_V->getData());
  // syevd_f('V',d_U,h_eigen);
  // Conditionnement
  Tcomp maxe = h_eigen->getData()[2];
  int cpt = 0;
  for (int i = 0; i < d_mat->getDims()[1]; i++) {
    if (h_eigen->getData()[i] < maxe / cond) {
      h_eigenvals_inv->getData()[i] = 0.;
      cpt++;
    } else
      h_eigenvals_inv->getData()[i] = 1. / h_eigen->getData()[i];
  }

  d_eigenvals_inv->host2device(h_eigenvals_inv->getData());

  carma_dgmm(cublas_handle, CUBLAS_SIDE_RIGHT, d_mat->getDims()[1],
             d_mat->getDims()[2], d_U->getData(), d_mat->getDims()[1],
             d_eigenvals_inv->getData(), 1, d_tmp->getData(),
             d_mat->getDims()[1]);
  carma_gemm<Tcomp>(cublas_handle, 'n', 't', d_mat->getDims()[1],
                    d_mat->getDims()[1], d_mat->getDims()[2], 1.0f,
                    d_tmp->getData(), d_mat->getDims()[1], d_U->getData(),
                    d_mat->getDims()[1], 0.0f, d_mat->getData(),
                    d_mat->getDims()[1]);

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
int sutra_controller_mv<Tcomp, Tout>::DDiago(carma_obj<Tcomp> *d_statcov,
                                             carma_obj<Tcomp> *d_geocov) {
  this->current_context->set_activeDevice(this->device, 1);
  const long dims_data[3] = {2, this->nactu(), this->nactu()};
  carma_obj<Tcomp> *d_M1 =
      new carma_obj<Tcomp>(this->current_context, dims_data);
  carma_obj<Tcomp> *d_tmp =
      new carma_obj<Tcomp>(this->current_context, dims_data);
  carma_obj<Tcomp> *d_tmp2 =
      new carma_obj<Tcomp>(this->current_context, dims_data);

  const long dims_data2[2] = {1, this->nactu()};
  carma_obj<Tcomp> *d_eigenvals =
      new carma_obj<Tcomp>(this->current_context, dims_data2);
  carma_obj<Tcomp> *d_eigenvals_sqrt =
      new carma_obj<Tcomp>(this->current_context, dims_data2);
  carma_obj<Tcomp> *d_eigenvals_inv =
      new carma_obj<Tcomp>(this->current_context, dims_data2);
  carma_host_obj<Tcomp> *h_eigenvals =
      new carma_host_obj<Tcomp>(dims_data2, MA_PAGELOCK);
  carma_host_obj<Tcomp> *h_eigenvals_inv =
      new carma_host_obj<Tcomp>(dims_data2, MA_PAGELOCK);
  carma_host_obj<Tcomp> *h_eigenvals_sqrt =
      new carma_host_obj<Tcomp>(dims_data2, MA_PAGELOCK);

  // 1. SVdec(geocov,U) --> Ut * geocov * U = D²
  carma_syevd<Tcomp>('V', d_geocov, h_eigenvals);

  d_eigenvals->host2device(h_eigenvals->getData());
  for (int i = 0; i < this->nactu(); i++) {
    h_eigenvals_sqrt->getData()[i] =
        sqrt(h_eigenvals->getData()[i]);  // D = sqrt(D²)
    h_eigenvals_inv->getData()[i] =
        1. / sqrt(h_eigenvals->getData()[i]);  // D⁻¹ = 1/sqrt(D²)
  }
  d_eigenvals_sqrt->host2device(h_eigenvals_sqrt->getData());
  d_eigenvals_inv->host2device(h_eigenvals_inv->getData());

  // 2. M⁻¹ = sqrt(eigenvals) * Ut : here, we have transpose(M⁻¹)

  carma_dgmm<Tcomp>(cublas_handle, CUBLAS_SIDE_RIGHT, this->nactu(),
                    this->nactu(), d_geocov->getData(), this->nactu(),
                    d_eigenvals_sqrt->getData(), 1, d_M1->getData(),
                    this->nactu());

  // 3. C' = M⁻¹ * statcov * M⁻¹t
  carma_gemm<Tcomp>(cublas_handle, 't', 'n', this->nactu(), this->nactu(),
                    this->nactu(), 1.0f, d_M1->getData(), this->nactu(),
                    d_statcov->getData(), this->nactu(), 0.0f, d_tmp->getData(),
                    this->nactu());

  carma_gemm<Tcomp>(cublas_handle, 'n', 'n', this->nactu(), this->nactu(),
                    this->nactu(), 1.0f, d_tmp->getData(), this->nactu(),
                    d_M1->getData(), this->nactu(), 0.0f, d_tmp2->getData(),
                    this->nactu());

  // 4. SVdec(C',A)
  carma_syevd<Tcomp>('V', d_tmp2, h_eigenvals);

  // 5. M = U * D⁻¹
  carma_dgmm<Tcomp>(cublas_handle, CUBLAS_SIDE_RIGHT, this->nactu(),
                    this->nactu(), d_geocov->getData(), this->nactu(),
                    d_eigenvals_inv->getData(), 1, d_tmp->getData(),
                    this->nactu());

  // 6. B = M * A;
  carma_gemm<Tcomp>(cublas_handle, 'n', 'n', this->nactu(), this->nactu(),
                    this->nactu(), 1.0f, d_tmp->getData(), this->nactu(),
                    d_tmp2->getData(), this->nactu(), 0.0f,
                    d_KLbasis->getData(), this->nactu());

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
  this->current_context->set_activeDevice(this->device, 1);
  this->d_covmat->host2device(covmat);
  return EXIT_SUCCESS;
}
template <typename Tcomp, typename Tout>
int sutra_controller_mv<Tcomp, Tout>::load_klbasis(Tcomp *klbasis) {
  this->current_context->set_activeDevice(this->device, 1);
  this->d_KLbasis->host2device(klbasis);
  return EXIT_SUCCESS;
}

// Florian features
template <typename Tcomp, typename Tout>
int sutra_controller_mv<Tcomp, Tout>::build_cmat(Tcomp cond) {
  this->current_context->set_activeDevice(this->device, 1);
  if (this->h_Cmmeigenvals != 0L) delete this->h_Cmmeigenvals;

  long Nactu = this->d_Cphim->getDims()[1];

  // (Cmm + Cn)⁻¹
  add_md(this->d_Cmm->getData(), this->d_Cmm->getData(),
         this->d_noisemat->getData(), this->nslope(),
         this->current_context->get_device(this->device));
  // invgen(this->d_Cmm,/*(Tcomp)(this->nslope()-this->nactu())*/200.0f,0);
  long dims_data1[2] = {1, this->nslope()};
  this->h_Cmmeigenvals = new carma_host_obj<Tcomp>(dims_data1, MA_PAGELOCK);

  invgen(this->d_Cmm, this->h_Cmmeigenvals, cond);

  // Cphim * (Cmm + Cn)⁻¹
  carma_gemm(cublas_handle, 'n', 'n', Nactu, this->nslope(), this->nslope(),
             1.0f, this->d_Cphim->getData(), Nactu, this->d_Cmm->getData(),
             this->nslope(), 0.0f, d_cmat->getData(), this->nactu());

  return EXIT_SUCCESS;
}
template <typename Tcomp, typename Tout>
int sutra_controller_mv<Tcomp, Tout>::filter_cmat(Tcomp cond) {
  this->current_context->set_activeDevice(this->device, 1);

  if (this->h_eigenvals != 0L) delete this->h_eigenvals;
  long Nactu = this->d_Cphim->getDims()[1];
  if (Nactu < this->nactu()) {
    long *dims_data = new long[3];
    dims_data[0] = 2;

    dims_data[1] = this->nslope();
    dims_data[2] = 2;
    carma_obj<Tcomp> *d_M =
        new carma_obj<Tcomp>(this->current_context, dims_data);
    dims_data[1] = Nactu;
    dims_data[2] = 2;
    carma_obj<Tcomp> *d_TT2ho =
        new carma_obj<Tcomp>(this->current_context, dims_data);
    dims_data[1] = 2;
    dims_data[2] = this->nslope();
    carma_obj<Tcomp> *d_M1 =
        new carma_obj<Tcomp>(this->current_context, dims_data);
    dims_data[2] = 2;
    carma_obj<Tcomp> *d_tmp3 =
        new carma_obj<Tcomp>(this->current_context, dims_data);

    // Imat decomposition TT
    dims_data[1] = Nactu;
    dims_data[2] = Nactu;
    carma_obj<Tcomp> *d_tmp2 =
        new carma_obj<Tcomp>(this->current_context, dims_data);

    // Dm⁻¹
    carma_gemm(cublas_handle, 't', 'n', Nactu, Nactu, this->nslope(), 1.0f,
               d_imat->getData(), this->nslope(), d_imat->getData(),
               this->nslope(), 0.0f, d_tmp2->getData(), Nactu);

    long dims_data1[2] = {1, Nactu};
    this->h_eigenvals = new carma_host_obj<Tcomp>(dims_data1, MA_PAGELOCK);
    invgen(d_tmp2, this->h_eigenvals, cond);

    dims_data[1] = Nactu;
    dims_data[2] = this->nslope();
    carma_obj<Tcomp> *d_Dm1 =
        new carma_obj<Tcomp>(this->current_context, dims_data);
    carma_gemm(cublas_handle, 'n', 't', Nactu, this->nslope(), Nactu, 1.0f,
               d_tmp2->getData(), Nactu, d_imat->getData(), this->nslope(),
               0.0f, d_Dm1->getData(), Nactu);

    delete d_tmp2;

    // TT2ho = Dm⁻¹ * Dtt
    carma_gemm(cublas_handle, 'n', 'n', Nactu, 2, this->nslope(), 1.0f,
               d_Dm1->getData(), Nactu,
               d_imat->getDataAt(this->nslope() * (Nactu)), this->nslope(),
               0.0f, d_TT2ho->getData(), Nactu);

    delete d_Dm1;

    // M = Dm * TT2ho
    carma_gemm(cublas_handle, 'n', 'n', this->nslope(), 2, Nactu, 1.0f,
               d_imat->getData(), this->nslope(), d_TT2ho->getData(), Nactu,
               0.0f, d_M->getData(), this->nslope());

    // M⁻¹
    carma_gemm(cublas_handle, 't', 'n', 2, 2, this->nslope(), 1.0f,
               d_M->getData(), this->nslope(), d_M->getData(), this->nslope(),
               0.0f, d_tmp3->getData(), 2);
    invgen(d_tmp3, 0.0f, 0);

    carma_gemm(cublas_handle, 'n', 't', 2, this->nslope(), 2, 1.0f,
               d_tmp3->getData(), 2, d_M->getData(), this->nslope(), 0.0f,
               d_M1->getData(), 2);

    // M*M⁻¹
    dims_data[1] = this->nslope();
    dims_data[2] = this->nslope();
    carma_obj<Tcomp> *d_Ftt =
        new carma_obj<Tcomp>(this->current_context, dims_data);
    carma_gemm(cublas_handle, 'n', 'n', this->nslope(), this->nslope(), 2, 1.0f,
               d_M->getData(), this->nslope(), d_M1->getData(), 2, 0.0f,
               d_Ftt->getData(), this->nslope());

    // TTcomp filter
    TT_filt(d_Ftt->getData(), this->nslope(),
            this->current_context->get_device(this->device));

    // cmat without TT
    dims_data[1] = Nactu;
    dims_data[2] = this->nslope();
    carma_obj<Tcomp> *d_cmat_tt =
        new carma_obj<Tcomp>(this->current_context, dims_data);

    carma_gemm(cublas_handle, 'n', 'n', Nactu, this->nslope(), this->nslope(),
               1.0f, d_cmat->getData(), this->nactu(), d_Ftt->getData(),
               this->nslope(), 0.0f, d_cmat_tt->getData(), Nactu);

    delete d_Ftt;

    // Fill CMAT
    fill_cmat(this->d_cmat->getData(), d_cmat_tt->getData(), d_M1->getData(),
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

  this->current_context->set_activeDevice(this->device, 1);
  if (strcmp(method, "inv") == 0) {
    //  R = (Dt*Cn⁻¹*D + Cphi⁻¹)⁻¹*Dt*Cn⁻¹

    carma_obj<Tcomp> *d_tmp;
    carma_obj<Tcomp> *d_tmp2;
    carma_obj<Tcomp> *d_tmp3;
    long *dims_data2 = new long[3];
    dims_data2[0] = 2;
    dims_data2[1] = this->nactu();
    dims_data2[2] = this->nactu();
    d_tmp = new carma_obj<Tcomp>(this->current_context, dims_data2);
    d_tmp2 = new carma_obj<Tcomp>(this->current_context, dims_data2);
    dims_data2[1] = this->nslope();
    dims_data2[2] = this->nactu();
    d_tmp3 = new carma_obj<Tcomp>(this->current_context, dims_data2);

    carma_dgmm(cublas_handle, CUBLAS_SIDE_LEFT, this->nslope(), this->nactu(),
               d_imat->getData(), this->nslope(), d_noisemat->getData(), 1,
               d_tmp3->getData(), this->nslope());
    carma_gemm(cublas_handle, 't', 'n', this->nactu(), this->nactu(),
               this->nslope(), one, d_tmp3->getData(), this->nslope(),
               d_imat->getData(), this->nslope(), zero, d_tmp2->getData(),
               this->nactu());
    delete d_tmp3;

    carma_geam(cublas_handle, 'n', 'n', this->nactu(), this->nactu(), one,
               d_tmp2->getData(), this->nactu(), one, d_covmat->getData(),
               this->nactu(), d_tmp->getData(), this->nactu());
    delete d_tmp2;

    carma_potri<Tcomp>(d_tmp);

    dims_data2[1] = this->nactu();
    dims_data2[2] = this->nslope();
    d_tmp2 = new carma_obj<Tcomp>(this->current_context, dims_data2);
    carma_gemm(cublas_handle, 'n', 't', this->nactu(), this->nslope(),
               this->nactu(), one, d_tmp->getData(), this->nactu(),
               d_imat->getData(), this->nslope(), zero, d_tmp2->getData(),
               this->nactu());
    delete d_tmp;

    carma_dgmm(cublas_handle, CUBLAS_SIDE_RIGHT, this->nactu(), this->nslope(),
               d_tmp2->getData(), this->nactu(), d_noisemat->getData(), 1,
               d_cmat->getData(), this->nactu());

    delete d_tmp2;
  }

  else if (strcmp(method, "n") == 0) {
    //  R = Cphi*Dt*(D*Cphi*Dt + Cn)⁻¹

    carma_obj<Tcomp> *d_tmp;
    carma_obj<Tcomp> *d_tmp2;
    carma_obj<Tcomp> *d_tmp3;
    carma_obj<Tcomp> *d_tmp4;
    long *dims_data2 = new long[3];
    dims_data2[0] = 2;
    dims_data2[1] = this->nslope();
    dims_data2[2] = this->nslope();
    d_tmp = new carma_obj<Tcomp>(this->current_context, dims_data2);
    d_tmp2 = new carma_obj<Tcomp>(this->current_context, dims_data2);
    d_tmp4 = new carma_obj<Tcomp>(this->current_context, dims_data2);
    //    carma_obj<Tcomp> *d_U;
    //    d_U = new carma_obj<Tcomp>(this->current_context, dims_data2);
    dims_data2[1] = this->nslope();
    dims_data2[2] = this->nactu();
    d_tmp3 = new carma_obj<Tcomp>(this->current_context, dims_data2);
    long *dims_data = new long[2];
    dims_data[0] = 1;
    dims_data[1] = this->nactu();
    // carma_host_obj<Tcomp> *h_eigenvals = new carma_host_obj<Tcomp>(dims_data,
    // MA_PAGELOCK); carma_obj<Tcomp> *d_eigenvals = new
    // carma_obj<Tcomp>(this->current_context, dims_data);

    carma_gemm(cublas_handle, 'n', 'n', this->nslope(), this->nactu(),
               this->nactu(), one, d_imat->getData(), this->nslope(),
               d_covmat->getData(), this->nactu(), zero, d_tmp3->getData(),
               this->nslope());
    carma_gemm(cublas_handle, 'n', 't', this->nslope(), this->nslope(),
               this->nactu(), one, d_tmp3->getData(), this->nslope(),
               d_imat->getData(), this->nslope(), zero, d_tmp2->getData(),
               this->nslope());
    delete d_tmp3;
    add_md(d_tmp4->getData(), d_tmp4->getData(), d_noisemat->getData(),
           this->nslope(), this->current_context->get_device(this->device));
    carma_geam(cublas_handle, 'n', 'n', this->nslope(), this->nslope(), one,
               d_tmp2->getData(), this->nslope(), one, d_tmp4->getData(),
               this->nslope(), d_tmp->getData(), this->nslope());
    delete d_tmp2;
    delete d_tmp4;

    carma_potri<Tcomp>(d_tmp);

    dims_data2[1] = this->nactu();
    dims_data2[2] = this->nslope();
    d_tmp2 = new carma_obj<Tcomp>(this->current_context, dims_data2);
    carma_gemm(cublas_handle, 't', 'n', this->nactu(), this->nslope(),
               this->nslope(), one, d_imat->getData(), this->nslope(),
               d_tmp->getData(), this->nslope(), zero, d_tmp2->getData(),
               this->nactu());
    delete d_tmp;
    carma_gemm(cublas_handle, 'n', 'n', this->nactu(), this->nslope(),
               this->nactu(), one, d_covmat->getData(), this->nactu(),
               d_tmp2->getData(), this->nactu(), zero, d_cmat->getData(),
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

  this->current_context->set_activeDevice(this->device, 1);
  if ((int)this->delay > 0) {
    for (int cc = 0; cc < this->delay; cc++)
      shift_buf(&((this->d_cenbuff->getData())[cc * this->nslope()]), 1,
                this->nslope(),
                this->current_context->get_device(this->device));

    carmaSafeCall(cudaMemcpy(
        &(this->d_cenbuff->getData()[(int)this->delay * this->nslope()]),
        this->d_centroids->getData(), sizeof(Tcomp) * this->nslope(),
        cudaMemcpyDeviceToDevice));

    carmaSafeCall(
        cudaMemcpy(this->d_centroids->getData(), this->d_cenbuff->getData(),
                   sizeof(Tcomp) * this->nslope(), cudaMemcpyDeviceToDevice));
  }

  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int sutra_controller_mv<Tcomp, Tout>::comp_com() {
  this->current_context->set_activeDevice(this->device, 1);

  // this->frame_delay();

  // POLC equations

  carma_geam<Tcomp>(
      cublas_handle, 'n', 'n', this->nactu(), 1, (Tcomp)(this->delay - 1),
      this->d_com1->getData(), this->nactu(), 1.0f - (this->delay - 1),
      this->d_com->getData(), this->nactu(), *d_compbuff, this->nactu());
  carma_gemv<Tcomp>(cublas_handle, 'n', this->nslope(), this->nactu(), 1.0f,
                    *d_imat, this->nslope(), *d_compbuff, 1, 0.0f, *d_compbuff2,
                    1);
  carma_geam<Tcomp>(cublas_handle, 'n', 'n', this->nslope(), 1, 1.0f,
                    *this->d_centroids, this->nslope(), -1.0f, *d_compbuff2,
                    this->nslope(), *d_olmeas, this->nslope());

  int nstreams = this->streams->get_nbStreams();
  if (nstreams > 1) {
    Tcomp alpha = -1.0f;
    Tcomp beta = 0.0f;

    for (int i = 0; i < nstreams; i++) {
      int istart1 =
          i * this->d_cmat->getDims(2) * this->d_cmat->getDims(1) / nstreams;
      int istart2 = i * this->d_cmat->getDims(1) / nstreams;

      cublasSetStream(cublas_handle, this->streams->get_stream(i));

      cublasOperation_t trans = carma_char2cublasOperation('n');

      carma_checkCublasStatus(
          cublasSgemv(cublas_handle, trans, this->d_cmat->getDims(1) / nstreams,
                      this->d_cmat->getDims(2), &alpha,
                      &((this->d_cmat->getData())[istart1]),
                      this->d_cmat->getDims(1) / nstreams, *d_olmeas, 1, &beta,
                      &((this->d_err->getData())[istart2]), 1));
    }

    this->streams->wait_all_streams();

  } else {
    // compute error
    this->d_err->gemv('n', -1.0f, this->d_cmat, this->d_cmat->getDims(1),
                      d_olmeas, 1, 0.0f, 1);  // POLC --> d_olmeas
  }
  /*
   mult_int(this->d_com->getData(), this->d_err->getData(),
   this->d_gain->getData(), this->gain, this->nactu(),
   this->current_context->get_device(device));*/

  carma_geam<Tcomp>(cublas_handle, 'n', 'n', this->nactu(), 1, this->gain,
                    *d_err, this->nactu(), 1.0f - this->gain, *this->d_com1,
                    this->nactu(), *this->d_com, this->nactu());

  return EXIT_SUCCESS;
}

template class sutra_controller_mv<float, float>;
template class sutra_controller_mv<float, uint16_t>;
