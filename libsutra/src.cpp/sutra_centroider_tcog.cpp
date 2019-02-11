#include <sutra_centroider_tcog.h>
#include <string>

template <class Tin, class T>
sutra_centroider_tcog<Tin, T>::sutra_centroider_tcog(carma_context *context,
                                                     sutra_wfs *wfs,
                                                     long nvalid, float offset,
                                                     float scale, int device)
    : sutra_centroider<Tin, T>(context, wfs, nvalid, offset, scale, device) {
  context->set_activeDevice(device, 1);

  this->nslopes = 2 * nvalid;
  this->threshold = 0;
  long dims_data2[2] = {1, this->nslopes};
  this->d_centroids_ref = new carma_obj<T>(this->current_context, dims_data2);
  this->d_centroids_ref->reset();
}

template <class Tin, class T>
sutra_centroider_tcog<Tin, T>::~sutra_centroider_tcog() {}

template <class Tin, class T>
string sutra_centroider_tcog<Tin, T>::get_type() {
  return "tcog";
}

template <class Tin, class T>
int sutra_centroider_tcog<Tin, T>::set_threshold(T threshold) {
  this->threshold = threshold;

  return EXIT_SUCCESS;
}

template <class Tin, class T>
int sutra_centroider_tcog<Tin, T>::get_cog(T *img, T *intensities, T *centroids,
                                           int nvalid, int npix, int ntot) {
  this->current_context->set_activeDevice(this->device, 1);

  get_centroids(ntot, (npix * npix), nvalid, npix, img, centroids,
                this->d_centroids_ref->getData(), this->d_validx->getData(),
                this->d_validy->getData(), intensities, this->threshold,
                this->scale, this->offset,
                this->current_context->get_device(this->device));

  // TODO: Implement get_cog_async
  // #ifndef USE_OLD
  //   // Modif Nono !!!
  //   // Now subap_reduce modify cube to apply the threshold on the cube so
  //   // get_centroids should apply threshold a second time
  //   subap_reduce_new(ntot, npix * npix, nvalid, cube, intensities,
  //                    this->threshold,
  //                    this->current_context->get_device(device));

  //   get_centroids(ntot, npix * npix, nvalid, npix, cube, centroids,
  //   intensities,
  //                 this->threshold, this->scale, this->offset,
  //                 this->current_context->get_device(device));
  // #else
  //   subap_reduce(ntot, npix * npix, nvalid, cube, intensities,
  //   this->threshold,
  //                this->current_context->get_device(device));

  //   get_centroids(ntot, npix * npix, nvalid, npix, cube, centroids,
  //   intensities,
  //                 this->threshold, this->scale, this->offset,
  //                 this->current_context->get_device(device));
  // #endif
  return EXIT_SUCCESS;
}

template <class Tin, class T>
int sutra_centroider_tcog<Tin, T>::get_cog(T *intensities, T *slopes,
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
int sutra_centroider_tcog<Tin, T>::get_cog() {
  if (this->wfs != nullptr)
    return this->get_cog(*(this->wfs->d_intensities), *(this->wfs->d_slopes),
                         true);
  DEBUG_TRACE("this->wfs was not initialized");
  return EXIT_FAILURE;
}

template class sutra_centroider_tcog<float, float>;
template class sutra_centroider_tcog<uint16_t, float>;
