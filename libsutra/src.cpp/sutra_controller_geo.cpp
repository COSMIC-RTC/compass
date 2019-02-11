#include <sutra_controller_geo.h>

template <typename T, typename Tout>
sutra_controller_geo<T, Tout>::sutra_controller_geo(carma_context *context,
                                                    long nactu, long Nphi,
                                                    float delay, sutra_dms *dms,
                                                    int *idx_dms, int ndm,
                                                    bool wfs_direction)
    : sutra_controller<T, Tout>(context, 0, 0, nactu, 0.0f, dms, idx_dms, ndm) {
  this->gain = 0.0f;
  this->Nphi = Nphi;

  //	long dims_data2[3];
  //	dims_data2[0] = 2;
  //	dims_data2[1] = nactu;
  //	dims_data2[2] = Nphi;
  // this->d_proj = new carma_obj<float>(this->current_context, dims_data2);
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
      this->d_cenbuff = new carma_obj<float>(this->current_context, dims_data2);
    }
    */
  this->Ntt = 0;

  long dims_data1[2];
  dims_data1[0] = 1;
  dims_data1[1] = nactu;
  this->d_gain = new carma_obj<T>(this->current_context, dims_data1);
  this->d_compfloat = new carma_obj<float>(this->current_context, dims_data1);
  this->d_compdouble = new carma_obj<double>(this->current_context, dims_data1);
  dims_data1[1] = Nphi;
  this->d_phi = new carma_obj<double>(this->current_context, dims_data1);
  this->d_indx_pup = new carma_obj<int>(this->current_context, dims_data1);
  if (wfs_direction)
    this->d_indx_mpup = new carma_obj<int>(this->current_context, dims_data1);
  else
    this->d_indx_mpup = 0L;
}

template <typename T, typename Tout>
sutra_controller_geo<T, Tout>::~sutra_controller_geo() {
  this->current_context->set_activeDevice(this->device, 1);
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
int sutra_controller_geo<T, Tout>::set_gain(T gain) {
  this->gain = gain;
  return EXIT_SUCCESS;
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
  this->d_geocov = new carma_obj<T>(this->current_context, dims_data);
  */
  this->d_geocov->host2device(Btt_pzt);
  this->d_geocovTT->host2device(Btt_TT);

  return EXIT_SUCCESS;
}
template <typename T, typename Tout>
int sutra_controller_geo<T, Tout>::init_proj(sutra_dms *dms, int *indx_dm,
                                             T *unitpervolt, int *indx_pup) {
  this->current_context->set_activeDevice(this->device, 1);
  long dims_data[3] = {2, this->Nphi, this->nactu()};
  carma_obj<T> d_IF(this->current_context, dims_data);
  dims_data[1] = this->nactu();
  carma_obj<T> d_tmp(this->current_context, dims_data);
  long dims_data1[2] = {1, this->Nphi * dms->ndm()};
  carma_obj<int> d_indx(this->current_context, dims_data1, indx_dm);

  this->d_indx_pup->host2device(indx_pup);

  // Get influence functions in d_IF
  int indx_start = 0;
  int ind = 0;
  vector<sutra_dm *>::iterator p;
  p = dms->d_dms.begin();
  while (p != dms->d_dms.end()) {
    sutra_dm *dm = *p;
    dm->get_IF<T>(d_IF.getDataAt(indx_start * this->Nphi),
                  d_indx.getDataAt(this->Nphi * ind), this->Nphi,
                  1.0f /*unitpervolt[ind]*/);
    indx_start += dm->nactus;
    ind++;
    p++;
  }

  // d_tmp = (transpose(d_IF) * d_IF)⁻¹
  carma_gemm<T>(this->cublas_handle(), 't', 'n', this->nactu(), this->nactu(),
                this->Nphi, 1.0f, d_IF.getData(), d_IF.getDims()[1],
                d_IF.getData(), d_IF.getDims()[1], 0.0f, d_tmp.getData(),
                d_tmp.getDims()[1]);
  carma_potri(&d_tmp);
  // invgen(d_tmp,1000.0f,1);

  // d_proj = d_tmp * transpose(d_IF)
  carma_gemm<T>(this->cublas_handle(), 'n', 't', this->nactu(), this->Nphi,
                this->nactu(), 1.0f, d_tmp.getData(), d_tmp.getDims()[1],
                d_IF.getData(), d_IF.getDims()[1], 0.0f,
                this->d_proj->getData(), this->d_proj->getDims()[1]);

  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int sutra_controller_geo<T, Tout>::init_proj_sparse(
    sutra_dms *dms, int *indx_dm, T *unitpervolt, int *indx_pup, int *indx_mpup,
    bool roket) {
  this->current_context->set_activeDevice(this->device, 1);
  vector<sutra_dm *>::iterator p;
  if (roket) {
    this->Ntt = 0;
    p = dms->d_dms.begin();
    while (p != dms->d_dms.end()) {
      sutra_dm *dm = *p;
      if (dm->type == "tt") this->Ntt += 1;
      p++;
    }
  }

  int Npzt = dms->ndm() - this->Ntt;
  carma_sparse_obj<double> *d_IFi[Npzt];
  long dims_data1[2] = {1, dms->ndm() * this->Nphi};
  carma_obj<int> d_indx(this->current_context, dims_data1, indx_dm);

  this->d_indx_pup->host2device(indx_pup);
  if (this->d_indx_mpup != 0L) this->d_indx_mpup->host2device(indx_mpup);

  // Get influence functions of the DM #ind in d_IFi
  int indx_start = 0;
  int ind = 0;
  int nnz = 0;
  int NNZ[Npzt];
  int Nact[Npzt];

  p = dms->d_dms.begin();
  while (p != dms->d_dms.end() - this->Ntt) {
    sutra_dm *dm = *p;
    dm->get_IF_sparse<double>(d_IFi[ind], d_indx.getDataAt(this->Nphi * ind),
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
  carma_obj<double> d_val(this->current_context, dims_data);
  carma_obj<int> d_col(this->current_context, dims_data);
  dims_data[1] = (this->nactu() - 2 * this->Ntt) + 1;
  carma_obj<int> d_row(this->current_context, dims_data);
  int cpt[Npzt];
  int nact = 0;
  cpt[0] = 0;
  p = dms->d_dms.begin();

  for (int i = 0; i < Npzt; i++) {
    sutra_dm *dm = *p;
    carmaSafeCall(cudaMemcpyAsync(d_val.getDataAt(cpt[i]), d_IFi[i]->d_data,
                                  sizeof(double) * d_IFi[i]->nz_elem,
                                  cudaMemcpyDeviceToDevice));
    carmaSafeCall(cudaMemcpyAsync(d_col.getDataAt(cpt[i]), d_IFi[i]->d_colind,
                                  sizeof(int) * d_IFi[i]->nz_elem,
                                  cudaMemcpyDeviceToDevice));
    if (i == 0)
      carmaSafeCall(cudaMemcpyAsync(d_row.getData(), d_IFi[i]->d_rowind,
                                    sizeof(int) * (dm->nactus + 1),
                                    cudaMemcpyDeviceToDevice));
    else
      carmaSafeCall(cudaMemcpyAsync(
          d_row.getDataAt(nact + 1), &(d_IFi[i]->d_rowind[1]),
          sizeof(int) * (dm->nactus), cudaMemcpyDeviceToDevice));
    cpt[i + 1] = cpt[i] + d_IFi[i]->nz_elem;
    nact += dm->nactus;
    p++;
    delete d_IFi[i];
  }

  if (Npzt > 1) {
    dims_data[1] = Npzt;
    carma_obj<int> d_NNZ(this->current_context, dims_data);
    carma_obj<int> d_nact(this->current_context, dims_data);
    d_NNZ.host2device(NNZ);
    d_nact.host2device(Nact);

    adjust_csr_index(d_row.getDataAt(1), d_NNZ.getData(), d_nact.getData(),
                     this->nactu() - 2 * this->Ntt, Nact[0],
                     this->current_context->get_device(this->device));
  }

  long dims_data2[3] = {2, (dms->nact_total() - 2 * this->Ntt), this->Nphi};
  this->d_IFsparse = new carma_sparse_obj<double>(
      this->current_context, dims_data2, d_val.getData(), d_col.getData(),
      d_row.getData(), nnz, false);

  // d_geocov = (transpose(d_IF) * d_IF)⁻¹
  carma_sparse_obj<double> *d_tmp =
      new carma_sparse_obj<double>(this->current_context);
  dims_data2[2] = (dms->nact_total() - 2 * this->Ntt);
  this->d_geocov = new carma_obj<T>(this->current_context, dims_data2);
  carma_obj<double> *d_tmp2 =
      new carma_obj<double>(this->current_context, this->d_geocov->getDims());

  carma_gemm<double>(cusparse_handle(), 'n', 't', this->d_IFsparse,
                     this->d_IFsparse, d_tmp);
  carma_csr2dense<double>(d_tmp, d_tmp2->getData());
  doubletofloat(d_tmp2->getData(), this->d_geocov->getData(),
                this->d_geocov->getNbElem(),
                this->current_context->get_device(this->device));

  this->d_geocov->scale(1.0f / this->Nphi, 1);
  carma_potri(d_geocov);
  // invgen(d_geocov,2.0f,0);
  delete d_tmp;
  delete d_tmp2;
  if (this->Ntt) {
    dims_data2[1] = this->Nphi;     // 2*this->Ntt;
    dims_data2[2] = 2 * this->Ntt;  // this->Nphi;
    this->d_TT = new carma_obj<T>(this->current_context, dims_data2);
    dims_data2[1] = 2 * this->Ntt;
    this->d_geocovTT = new carma_obj<T>(this->current_context, dims_data2);
    dims_data1[1] = this->Nphi;
    this->d_phif = new carma_obj<T>(this->current_context, dims_data1);

    p = dms->d_dms.begin();
    ind = 0;
    int ind2 = 0;
    while (p != dms->d_dms.end()) {
      sutra_dm *dm = *p;
      if (dm->type == "tt") {
        // dm->get_IF(this->d_TT->getDataAt(ind*Nphi),
        // d_indx.getDataAt(this->Nphi * ind2), this->Nphi, 1.0f);
        for (int i = 0; i < dm->nactus; i++) {
          dm->comp_oneactu(i, 1.0f);

          getIF<T>(this->d_TT->getDataAt(ind * this->Nphi),
                   dm->d_shape->d_screen->getData(),
                   d_indx.getDataAt(ind2 * this->Nphi), this->Nphi, 0,
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
               this->Nphi, 1.0f / this->Nphi, this->d_TT->getData(), this->Nphi,
               this->d_TT->getData(), this->Nphi, 0.0f,
               this->d_geocovTT->getData(), 2 * this->Ntt);

    T *tmp;
    tmp = (T *)malloc(this->d_geocovTT->getNbElem() * sizeof(T));
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
int sutra_controller_geo<T, Tout>::comp_dphi(sutra_source *target,
                                             bool wfs_direction) {
  this->current_context->set_activeDevice(this->device, 1);
  // Get the target phase in the pupil
  if (wfs_direction && this->d_indx_mpup == 0L) {
    DEBUG_TRACE("controller geo has not been initialized for wfs direction");
    return EXIT_FAILURE;
  }
  if (wfs_direction)
    get_pupphase(this->d_phi->getData(), target->d_phase->d_screen->getData(),
                 this->d_indx_mpup->getData(), this->Nphi,
                 this->current_context->get_device(this->device));
  else
    get_pupphase(this->d_phi->getData(), target->d_phase->d_screen->getData(),
                 this->d_indx_pup->getData(), this->Nphi,
                 this->current_context->get_device(this->device));

  remove_avg(this->d_phi->getData(), this->d_phi->getNbElem(),
             this->current_context->get_device(this->device));

  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int sutra_controller_geo<T, Tout>::comp_com() {
  // Project the phase on the actuators
  /*
  //Dense version
  carma_gemv(this->cublas_handle(),'n', this->nactu(), this->Nphi, -1.0f,
  this->d_proj->getData(),this->d_proj->getDims()[1], this->d_phi->getData(),1,
  0.0f, this->d_com->getData(),1);
  */
  this->current_context->set_activeDevice(this->device, 1);

  // Sparse version
  carma_gemv(cusparse_handle(), 'n', 1.0 / this->Nphi, this->d_IFsparse,
             this->d_phi->getData(), 0.0, this->d_compdouble->getData());
  doubletofloat(this->d_compdouble->getData(), this->d_compfloat->getData(),
                this->nactu(), this->current_context->get_device(this->device));
  // If we are in error budget case, d_geocov is Btt*Btt.transpose
  carma_gemv(this->cublas_handle(), 'n', this->nactu() - 2 * this->Ntt,
             this->nactu() - 2 * this->Ntt, -1.0f, this->d_geocov->getData(),
             this->d_geocov->getDims()[1], this->d_compfloat->getData(), 1,
             0.0f, this->d_com->getData(), 1);
  if (this->Ntt) {
    doubletofloat(this->d_phi->getData(), this->d_phif->getData(), this->Nphi,
                  this->current_context->get_device(this->device));
    carma_gemv(this->cublas_handle(), 't', this->d_TT->getDims(1),
               this->d_TT->getDims(2), 1.0f / this->Nphi, this->d_TT->getData(),
               this->d_TT->getDims(1), this->d_phif->getData(), 1, 0.0f,
               this->d_compfloat->getData(), 1);
    carma_gemv(this->cublas_handle(), 'n', 2 * this->Ntt, 2 * this->Ntt, -1.0f,
               this->d_geocovTT->getData(), 2 * this->Ntt,
               this->d_compfloat->getData(), 1, 0.0f,
               this->d_com->getDataAt(this->d_IFsparse->getDims(1)), 1);
  }

  return EXIT_SUCCESS;
}

template class sutra_controller_geo<float, float>;
template class sutra_controller_geo<float, uint16_t>;
