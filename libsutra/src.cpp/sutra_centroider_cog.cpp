#include <sutra_centroider_cog.h>
#include <string>

template <typename T>
sutra_centroider_cog<T>::sutra_centroider_cog(carma_context *context,
                                              sutra_wfs *wfs, long nvalid,
                                              T offset, T scale, int device)
    : sutra_centroider<T>(context, wfs, nvalid, offset, scale, device) {
  this->nslopes = 2 * nvalid;
  long dims_data2[2] = {1, this->nslopes};
  this->d_centroids_ref = new carma_obj<T>(this->current_context, dims_data2);
  this->d_centroids_ref->reset();
}

template <typename T>
sutra_centroider_cog<T>::~sutra_centroider_cog() {}

template <typename T>
string sutra_centroider_cog<T>::get_type() {
  return "cog";
}

template <typename T>
int sutra_centroider_cog<T>::get_cog(T *img, T *intensities, T *centroids,
                                     int nvalid, int npix, int ntot) {
  this->current_context->set_activeDevice(this->device, 1);

  // subap_reduce(ntot, (npix * npix), nvalid, img, ref,
  //              this->current_context->get_device(device));
  get_centroids(ntot, (npix * npix), nvalid, npix, img, centroids,
                this->d_centroids_ref->getData(), this->d_validx->getData(),
                this->d_validy->getData(), intensities, this->scale,
                this->offset, this->current_context->get_device(this->device));

  return EXIT_SUCCESS;
}

template <typename T>
int sutra_centroider_cog<T>::get_cog(T *intensities, T *slopes, bool noise) {
  if (this->wfs != nullptr) {
    if (noise || this->wfs->roket == false) {
      return this->get_cog(*(this->wfs->d_binimg), intensities, slopes,
                           this->wfs->nvalid_tot, this->wfs->npix,
                           this->wfs->d_binimg->getDims()[1]);
    } else {
      return this->get_cog(*(this->wfs->d_binimg_notnoisy), intensities, slopes,
                           this->wfs->nvalid_tot, this->wfs->npix,
                           this->wfs->d_binimg->getDims()[1]);
    }
  }
  DEBUG_TRACE("this->wfs was not initialized");
  return EXIT_SUCCESS;
}

template <typename T>
int sutra_centroider_cog<T>::get_cog() {
  if (this->wfs != nullptr)
    return this->get_cog(*(this->wfs->d_intensities), *(this->wfs->d_slopes),
                         true);

  DEBUG_TRACE("this->wfs was not initialized");
  return EXIT_SUCCESS;
}

template class sutra_centroider_cog<float>;

#ifdef CAN_DO_HALF
template <>
int sutra_centroider_cog<half>::get_cog(half *intensities, half *slopes,
                                        bool noise) {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_FAILURE;
}

template <>
int sutra_centroider_cog<half>::get_cog() {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_FAILURE;
}

template class sutra_centroider_cog<half>;
#endif
