#include <sutra_centroider_wcog.h>
#include <string>

template <class Tin, class T>
sutra_centroider_wcog<Tin, T>::sutra_centroider_wcog(carma_context *context,
                                                     sutra_wfs *wfs,
                                                     long nvalid, float offset,
                                                     float scale, bool filter_TT,
						     int device)
  : sutra_centroider<Tin, T>(context, wfs, nvalid, offset, scale, filter_TT, device) {
  context->set_activeDevice(device, 1);

  this->nslopes = 2 * nvalid;
  if (wfs != nullptr)
    this->npix = this->wfs->npix;
  else
    this->npix = 0;
  this->d_weights = 0L;
  long dims_data2[2] = {1, this->nslopes};
  this->d_centroids_ref = new carma_obj<T>(this->current_context, dims_data2);
  this->d_centroids_ref->reset();

  if(filter_TT){
    this->init_TT_filter();
  }
  
}

template <class Tin, class T>
sutra_centroider_wcog<Tin, T>::~sutra_centroider_wcog() {}

template <class Tin, class T>
string sutra_centroider_wcog<Tin, T>::get_type() {
  return "wcog";
}

template <class Tin, class T>
int sutra_centroider_wcog<Tin, T>::init_weights() {
  this->current_context->set_activeDevice(this->device, 1);
  if (this->d_weights != 0L) delete this->d_weights;

  long *dims_data3 = new long[4];
  dims_data3[0] = 3;
  dims_data3[1] = this->npix;
  dims_data3[2] = this->npix;
  dims_data3[3] = this->nvalid;

  this->current_context->set_activeDevice(this->device, 1);
  this->d_weights = new carma_obj<float>(this->current_context, dims_data3);

  delete[] dims_data3;

  return EXIT_SUCCESS;
}

template <class Tin, class T>
int sutra_centroider_wcog<Tin, T>::load_weights(float *weights, int ndim) {
  if (ndim == 3)
    this->d_weights->host2device(weights);
  else {
    // weights is a 2d array
    // same weight for each subap
    float *tmp;  ///< Input data
    carmaSafeCall(
        cudaMalloc((void **)&tmp, sizeof(float) * this->npix * this->npix));
    carmaSafeCall(cudaMemcpy(tmp, weights,
                             sizeof(float) * this->npix * this->npix,
                             cudaMemcpyHostToDevice));
    fillweights<float>(*(this->d_weights), tmp, this->npix,
                       this->d_weights->getNbElem(),
                       this->current_context->get_device(this->device));
    carmaSafeCall(cudaFree(tmp));
  }

  return EXIT_SUCCESS;
}

template <class Tin, class T>
int sutra_centroider_wcog<Tin, T>::set_npix(int npix) {
  this->npix = npix;
  return EXIT_SUCCESS;
}
template <class Tin, class T>
int sutra_centroider_wcog<Tin, T>::get_cog(float *img, float *intensities,
                                           T *centroids, int nvalid, int npix,
                                           int ntot) {
  // wcog
  // TODO: Implement sutra_centroider_wcog<Tin, T>::get_cog_async
  // subap_reduce<T>(ntot, npix * npix, nvalid, cube, intensities,
  //                     *(this->d_weights),
  //                     this->current_context->get_device(device));

  // get_centroids<T>(ntot, npix * npix, nvalid, npix, cube, centroids,
  //                      intensities, *(this->d_weights), this->scale,
  //                      this->offset,
  //                      this->current_context->get_device(device));

  get_centroids(ntot, (npix * npix), nvalid, npix, img, centroids,
                this->d_centroids_ref->getData(), this->d_validx->getData(),
                this->d_validy->getData(), intensities,
                this->d_weights->getData(), this->scale, this->offset,
                this->current_context->get_device(this->device));

  if (this->filter_TT) {
    this->apply_TT_filter(centroids);
 }

 return EXIT_SUCCESS;
}

template <class Tin, class T>
int sutra_centroider_wcog<Tin, T>::get_cog(float *intensities, T *slopes,
                                           bool noise) {
  if (this->wfs != nullptr) {
    if (noise || this->wfs->roket == false)
      return this->get_cog(*(this->wfs->d_binimg), intensities, slopes,
                           this->wfs->nvalid, this->wfs->npix,
                           this->wfs->d_binimg->getDims()[1]);
    else
      return this->get_cog(*(this->wfs->d_binimg_notnoisy), intensities, slopes,
                           this->wfs->nvalid, this->wfs->npix,
                           this->wfs->d_binimg->getDims()[1]);
  }
  DEBUG_TRACE("this->wfs was not initialized");
  return EXIT_FAILURE;
}

template <class Tin, class T>
int sutra_centroider_wcog<Tin, T>::get_cog() {
  if (this->wfs != nullptr)
    return this->get_cog(*(this->wfs->d_intensities), *(this->wfs->d_slopes),
                         true);
  DEBUG_TRACE("this->wfs was not initialized");
  return EXIT_FAILURE;
}

template class sutra_centroider_wcog<float, float>;
template class sutra_centroider_wcog<uint16_t, float>;
