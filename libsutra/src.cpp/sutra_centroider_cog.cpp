#include <sutra_centroider_cog.h>
#include <string>

template <class Tin, class T>
sutra_centroider_cog<Tin, T>::sutra_centroider_cog(carma_context *context,
                                                   sutra_wfs *wfs, long nvalid,
                                                   float offset, float scale,
                                                   bool filter_TT, int device)
  : sutra_centroider<Tin, T>(context, wfs, nvalid, offset, scale, filter_TT, device) {
  this->nslopes = 2 * nvalid;
  long dims_data2[2] = {1, this->nslopes};
  this->d_centroids_ref = new carma_obj<T>(this->current_context, dims_data2);
  this->d_centroids_ref->reset();
  
  long dims_data[2] = {1, this->nvalid};
  if (this->filter_TT == true) {
    dims_data[1] = 2;
    this->d_TT_slopes = new carma_obj<T>(this->current_context, dims_data);
    dims_data[1] = this->nslopes;
    this->d_ref_Tip = new carma_obj<T>(this->current_context, dims_data);
    this->d_ref_Tilt = new carma_obj<T>(this->current_context, dims_data);
  }
}

template <class Tin, class T>
sutra_centroider_cog<Tin, T>::~sutra_centroider_cog() {}

template <class Tin, class T>
string sutra_centroider_cog<Tin, T>::get_type() {
  return "cog";
}

template <class Tin, class T>
int sutra_centroider_cog<Tin, T>::get_cog(float *img, float *intensities,
                                          T *centroids, int nvalid, int npix,
                                          int ntot) {
  this->current_context->set_activeDevice(this->device, 1);

  // subap_reduce(ntot, (npix * npix), nvalid, img, ref,
  //              this->current_context->get_device(device));
  get_centroids(ntot, (npix * npix), nvalid, npix, img, centroids,
                this->d_centroids_ref->getData(), this->d_validx->getData(),
                this->d_validy->getData(), intensities, this->scale,
                this->offset, this->current_context->get_device(this->device));
  
  if (this->filter_TT == true) {
    this->wfs->d_slopes->copyFrom(centroids, this->nslopes);
    
    T tip = this->wfs->d_slopes->dot(this->d_ref_Tip,1,1);
    T tilt = this->wfs->d_slopes->dot(this->d_ref_Tilt,1,1);

    this->wfs->d_slopes->axpy(T(-1 * tip), this->d_ref_Tip, 1, 1);
    this->wfs->d_slopes->axpy(T(-1 * tilt), this->d_ref_Tilt, 1, 1);
    
    this->wfs->d_slopes->copyInto(centroids, this->nslopes);

 }
  
  
  return EXIT_SUCCESS;
}

template <class Tin, class T>
int sutra_centroider_cog<Tin, T>::get_cog(float *intensities, T *slopes,
                                          bool noise) {
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

template <class Tin, class T>
int sutra_centroider_cog<Tin, T>::get_cog() {
  if (this->wfs != nullptr)
    return this->get_cog(*(this->wfs->d_intensities), *(this->wfs->d_slopes),
                         true);

  DEBUG_TRACE("this->wfs was not initialized");
  return EXIT_SUCCESS;
}

template class sutra_centroider_cog<float, float>;
template class sutra_centroider_cog<uint16_t, float>;

#ifdef CAN_DO_HALF
template <>
int sutra_centroider_cog<float, half>::get_cog(float *intensities, half *slopes,
                                               bool noise) {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_FAILURE;
}

template <>
int sutra_centroider_cog<uint16_t, half>::get_cog(float *intensities,
                                                  half *slopes, bool noise) {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_FAILURE;
}

template <>
int sutra_centroider_cog<float, half>::get_cog() {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_FAILURE;
}

template <>
int sutra_centroider_cog<uint16_t, half>::get_cog() {
  DEBUG_TRACE("Not implemented for half precision");
  return EXIT_FAILURE;
}

template class sutra_centroider_cog<float, half>;
template class sutra_centroider_cog<uint16_t, half>;
#endif
