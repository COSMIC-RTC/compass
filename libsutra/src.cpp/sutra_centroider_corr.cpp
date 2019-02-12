#include <sutra_centroider_corr.h>
#include <string>

template <class Tin, class T>
sutra_centroider_corr<Tin, T>::sutra_centroider_corr(carma_context *context,
                                                     sutra_wfs *wfs,
                                                     long nvalid, float offset,
                                                     float scale, int device)
    : sutra_centroider<Tin, T>(context, wfs, nvalid, offset, scale, device) {
  context->set_activeDevice(device, 1);

  this->nslopes = 2 * nvalid;
  long dims_data2[2] = {1, this->nslopes};
  this->d_centroids_ref = new carma_obj<T>(this->current_context, dims_data2);
  this->d_centroids_ref->reset();

  this->d_corrfnct = 0L;
  this->d_corrspot = 0L;
  this->d_corrnorm = 0L;
  this->d_corrmax = 0L;
  this->d_corr = 0L;
  this->d_interpmat = 0L;

  if (wfs != nullptr) {
    this->npix = wfs->npix;
    this->nxsub = wfs->nxsub;
  } else {
    this->npix = 0;
    this->nxsub = 0;
  }
  long dims_data3[4] = {3, this->npix, this->npix, this->nvalid};
  this->d_bincube = new carma_obj<T>(this->current_context, dims_data3);

  this->interp_sizex = 0;
  this->interp_sizey = 0;
}

template <class Tin, class T>
sutra_centroider_corr<Tin, T>::~sutra_centroider_corr() {
  if (this->d_corrfnct != 0L) delete this->d_corrfnct;
  if (this->d_corrspot != 0L) delete this->d_corrspot;
  if (this->d_corrnorm != 0L) delete this->d_corrnorm;
  if (this->d_corrmax != 0L) delete this->d_corrmax;
  if (this->d_corr != 0L) delete this->d_corr;
  if (this->d_interpmat != 0L) delete this->d_interpmat;
}

template <class Tin, class T>
string sutra_centroider_corr<Tin, T>::get_type() {
  return "corr";
}

template <class Tin, class T>
int sutra_centroider_corr<Tin, T>::set_npix(int npix) {
  this->npix = npix;

  return EXIT_SUCCESS;
}

template <class Tin, class T>
int sutra_centroider_corr<Tin, T>::fill_bincube(T *img) {
  this->current_context->set_activeDevice(this->device, 1);

  fillbincube(img, this->d_bincube->getData(), this->npix, this->nvalid,
              this->nxsub, this->d_validx->getData(), this->d_validy->getData(),
              this->current_context->get_device(this->device));
  return EXIT_SUCCESS;
}

template <class Tin, class T>
int sutra_centroider_corr<Tin, T>::init_corr(int isizex, int isizey,
                                             T *interpmat) {
  this->current_context->set_activeDevice(this->device, 1);
  if (this->d_corrfnct != 0L) delete this->d_corrfnct;
  if (this->d_corrspot != 0L) delete this->d_corrspot;
  if (this->d_corrnorm != 0L) delete this->d_corrnorm;
  if (this->d_corrmax != 0L) delete this->d_corrmax;
  if (this->d_corr != 0L) delete this->d_corr;
  if (this->d_interpmat != 0L) delete this->d_interpmat;

  long *dims_data3 = new long[4];
  dims_data3[0] = 3;
  dims_data3[1] = 2 * this->npix;
  dims_data3[2] = 2 * this->npix;
  dims_data3[3] = this->nvalid;

  this->d_corrfnct =
      new carma_obj<cuFloatComplex>(this->current_context, dims_data3);
  this->d_corrspot =
      new carma_obj<cuFloatComplex>(this->current_context, dims_data3);

  int mdims[2];
  mdims[0] = (int)dims_data3[1];
  mdims[1] = (int)dims_data3[2];
  cufftHandle *plan = this->d_corrfnct->getPlan();  ///< FFT plan
  carmafftSafeCall(cufftPlanMany(plan, 2, mdims, NULL, 1, 0, NULL, 1, 0,
                                 CUFFT_C2C, (int)dims_data3[3]));

  dims_data3[1] = 2 * this->npix - 1;
  dims_data3[2] = 2 * this->npix - 1;
  this->d_corr = new carma_obj<T>(this->current_context, dims_data3);

  long *dims_data2 = new long[3];
  dims_data2[0] = 2;
  dims_data2[1] = 2 * this->npix - 1;
  dims_data2[2] = 2 * this->npix - 1;
  this->d_corrnorm = new carma_obj<T>(this->current_context, dims_data2);

  long *dims_data1 = new long[2];
  dims_data1[0] = 1;
  dims_data1[1] = this->nvalid;
  this->d_corrmax = new carma_obj<int>(this->current_context, dims_data1);

  this->interp_sizex = isizex;
  this->interp_sizey = isizey;
  dims_data2[1] = isizex * isizey;
  dims_data2[2] = 6;
  this->d_interpmat = new carma_obj<T>(this->current_context, dims_data2);
  this->d_interpmat->host2device(interpmat);

  delete[] dims_data1;
  delete[] dims_data2;
  delete[] dims_data3;

  return EXIT_SUCCESS;
}

template <class Tin, class T>
int sutra_centroider_corr<Tin, T>::load_corr(T *corr, T *corr_norm, int ndim) {
  this->current_context->set_activeDevice(this->device, 1);
  int nval = (ndim == 3) ? 1 : this->nvalid;

  this->d_corrnorm->host2device(corr_norm);

  T *tmp;  ///< Input data
  cudaError err;

  if (ndim == 3) {
    carmaSafeCall(err =
                      cudaMalloc((void **)&tmp, sizeof(T) * this->npix *
                                                    this->npix * this->nvalid));
    carmaSafeCall(
        err = cudaMemcpy(tmp, corr,
                         sizeof(T) * this->npix * this->npix * this->nvalid,
                         cudaMemcpyHostToDevice));
  } else {
    carmaSafeCall(
        err = cudaMalloc((void **)&tmp, sizeof(T) * this->npix * this->npix));
    carmaSafeCall(err =
                      cudaMemcpy(tmp, corr, sizeof(T) * this->npix * this->npix,
                                 cudaMemcpyHostToDevice));
  }

  fillcorr<cuFloatComplex, T>(*(this->d_corrfnct), tmp, this->npix,
                              2 * this->npix,
                              this->npix * this->npix * this->nvalid, nval,
                              this->current_context->get_device(this->device));

  carmaSafeCall(err = cudaFree(tmp));

  carma_fft<cuFloatComplex, cuFloatComplex>(*(this->d_corrfnct),
                                            *(this->d_corrfnct), 1,
                                            *this->d_corrfnct->getPlan());

  return EXIT_SUCCESS;
}

template <class Tin, class T>
int sutra_centroider_corr<Tin, T>::get_cog(float *img, T *intensities,
                                           T *centroids, int nvalid, int npix,
                                           int ntot) {
  this->current_context->set_activeDevice(this->device, 1);
  cudaError err;

  // set corrspot to 0
  carmaSafeCall(
      err = cudaMemset(*(this->d_corrspot), 0,
                       sizeof(cuFloatComplex) * this->d_corrspot->getNbElem()));

  if (this->nxsub != 0)
    fill_bincube(img);
  else
    std::cerr << "nxsub not set for this centroider" << std::endl;
  // correlation algorithm
  fillcorr<cuFloatComplex, T>(*(this->d_corrspot), this->d_bincube->getData(),
                              this->npix, 2 * this->npix,
                              this->npix * this->npix * this->nvalid, 1,
                              this->current_context->get_device(this->device));

  carma_fft<cuFloatComplex, cuFloatComplex>(*(this->d_corrspot),
                                            *(this->d_corrspot), 1,
                                            *this->d_corrfnct->getPlan());

  correl<cuFloatComplex>(*(this->d_corrspot), *(this->d_corrfnct),
                         this->d_corrfnct->getNbElem(),
                         this->current_context->get_device(this->device));
  // after this d_corrspot contains the fft of the correl function

  carma_fft<cuFloatComplex, cuFloatComplex>(*(this->d_corrspot),
                                            *(this->d_corrspot), -1,
                                            *this->d_corrfnct->getPlan());

  // size is 2 x npix so it is even ...
  roll2real<cuFloatComplex, T>(
      *(this->d_corr), *(this->d_corrspot), 2 * this->npix,
      (2 * this->npix) * (2 * this->npix), this->d_corrspot->getNbElem(),
      this->current_context->get_device(this->device));
  // here need to normalize
  corr_norm<T>(*(this->d_corr), *(this->d_corrnorm),
               this->d_corrnorm->getNbElem(), this->d_corr->getNbElem(),
               this->current_context->get_device(this->device));

  // need to find max for each subap
  // if the corr array for one subap is greater than 20x20
  // it won't fit in shared mem
  // so we window around the center of the array,
  // the max is expected to be found inside the npix x npix central part anyway
  int nbmax = (2 * this->npix - 1 > 20) ? this->npix : 2 * this->npix - 1;
  int xoff = this->d_corr->getDims(1) / 2 - nbmax / 2;
  int yoff = this->d_corr->getDims(2) / 2 - nbmax / 2;

  subap_sortmaxi<T>(nbmax * nbmax, this->nvalid, *(this->d_corr),
                    *(this->d_corrmax), 1, xoff, yoff, nbmax,
                    this->d_corr->getDims(1));

  // do parabolic interpolation
  subap_pinterp<T>(this->interp_sizex * this->interp_sizey, this->nvalid,
                   *(this->d_corr), *(this->d_corrmax), centroids,
                   *(this->d_interpmat), this->interp_sizex, this->interp_sizey,
                   this->nvalid, 2 * this->npix - 1, this->scale, this->offset);

  carma_axpy<T>(this->current_context->get_cublasHandle(), this->nslopes, -1.0f,
                this->d_centroids_ref->getData(), 1, centroids, 1);

  return EXIT_SUCCESS;
}

template <class Tin, class T>
int sutra_centroider_corr<Tin, T>::get_cog(T *intensities, T *slopes,
                                           bool noise) {
  if (this->wfs != nullptr) {
    return this->get_cog(*this->wfs->d_binimg, intensities, slopes,
                         this->wfs->nvalid, this->wfs->npix,
                         this->wfs->d_bincube->getNbElem());
  }

  DEBUG_TRACE("this->wfs was not initialized");
  return EXIT_FAILURE;
}

template <class Tin, class T>
int sutra_centroider_corr<Tin, T>::get_cog() {
  if (this->wfs != nullptr)
    return this->get_cog(*(this->wfs->d_intensities), *(this->wfs->d_slopes),
                         true);

  DEBUG_TRACE("this->wfs was not initialized");
  return EXIT_FAILURE;
}

template class sutra_centroider_corr<float, float>;
template class sutra_centroider_corr<uint16_t, float>;
