#include <sutra_centroider_bpcog.h>

template <class Tin, class T>
sutra_centroider_bpcog<Tin, T>::sutra_centroider_bpcog(
    carma_context *context, sutra_wfs *wfs, long nvalid, float offset,
    float scale, int device, int nmax)
    : sutra_centroider<Tin, T>(context, wfs, nvalid, offset, scale, device) {
  this->nslopes = 2 * nvalid;
  this->nmax = nmax;
  long dims_data2[2] = {1, this->nslopes};
  this->d_centroids_ref = new carma_obj<T>(this->current_context, dims_data2);
  this->d_centroids_ref->reset();

  long dims_data[3];
  dims_data[0] = 2;
  dims_data[1] = nmax;
  dims_data[2] = nvalid;

  this->d_bpix = new carma_obj<T>(this->current_context, dims_data);
  this->d_bpind = new carma_obj<uint>(this->current_context, dims_data);
}

template <class Tin, class T>
sutra_centroider_bpcog<Tin, T>::~sutra_centroider_bpcog() {
  delete this->d_bpix;
  delete this->d_bpind;
}

template <class Tin, class T>
string sutra_centroider_bpcog<Tin, T>::get_type() {
  return "bpcog";
}

template <class Tin, class T>
int sutra_centroider_bpcog<Tin, T>::set_nmax(int nmax) {
  this->current_context->set_activeDevice(this->device, 1);
  this->nmax = nmax;
  delete this->d_bpix;
  delete this->d_bpind;

  long dims_data[3];
  dims_data[0] = 2;
  dims_data[1] = nmax;
  dims_data[2] = this->nvalid;

  this->d_bpix = new carma_obj<T>(this->current_context, dims_data);
  this->d_bpind = new carma_obj<uint>(this->current_context, dims_data);

  return EXIT_SUCCESS;
}

template <class Tin, class T>
int sutra_centroider_bpcog<Tin, T>::get_cog(float *img, float *intensities,
                                            T *centroids, int nvalid, int npix,
                                            int ntot) {
  this->current_context->set_activeDevice(this->device, 1);

  get_centroids(ntot, (npix * npix), nvalid, npix, img, centroids,
                this->d_centroids_ref->getData(), this->d_validx->getData(),
                this->d_validy->getData(), intensities, this->nmax, this->scale,
                this->offset, this->current_context->get_device(this->device));

  // brightest pixels cog
  // subap_sortmax<T>(npix * npix, nvalid, cube, this->d_bpix->getData(),
  //                      this->d_bpind->getData(), this->nmax,
  //                      current_context->get_device(device));

  // subap_bpcentro<T>(this->nmax, nvalid, npix, this->d_bpix->getData(),
  //                       this->d_bpind->getData(), centroids, this->scale,
  //                       this->offset);

  return EXIT_SUCCESS;
}

template <class Tin, class T>
int sutra_centroider_bpcog<Tin, T>::get_cog(float *intensities, T *slopes,
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
int sutra_centroider_bpcog<Tin, T>::get_cog() {
  if (this->wfs != nullptr)
    return this->get_cog(*(this->wfs->d_intensities), *(this->wfs->d_slopes),
                         true);

  DEBUG_TRACE("this->wfs was not initialized");
  return EXIT_FAILURE;
}

template class sutra_centroider_bpcog<float, float>;
template class sutra_centroider_bpcog<uint16_t, float>;

#ifdef CAN_DO_HALF
template <>
int sutra_centroider_bpcog<float, half>::get_cog(float *intensities,
                                                 half *slopes, bool noise) {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_FAILURE;
}

template <>
int sutra_centroider_bpcog<uint16_t, half>::get_cog(float *intensities,
                                                    half *slopes, bool noise) {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_FAILURE;
}

template <>
int sutra_centroider_bpcog<float, half>::get_cog() {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_FAILURE;
}

template <>
int sutra_centroider_bpcog<uint16_t, half>::get_cog() {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_FAILURE;
}

template class sutra_centroider_bpcog<float, half>;
template class sutra_centroider_bpcog<uint16_t, half>;
#endif
