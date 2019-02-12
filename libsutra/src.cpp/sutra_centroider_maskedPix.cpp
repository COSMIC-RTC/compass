#include <sutra_centroider_maskedPix.h>
#include <string>

template <class Tin, class T>
sutra_centroider_maskedPix<Tin, T>::sutra_centroider_maskedPix(
    carma_context *context, sutra_wfs *wfs, long nvalid, long npupils,
    float offset, float scale, int device)
    : sutra_centroider<Tin, T>(context, wfs, nvalid, offset, scale, device) {
  context->set_activeDevice(device, 1);

  if (wfs != nullptr)
    this->nslopes = wfs->d_validsubsx->getNbElem();
  else
    this->nslopes = nvalid * npupils;
  long dims_data[2] = {1, this->nslopes};
  this->d_intensities = new carma_obj<T>(context, dims_data);
  this->d_intensities->init_reduceCub();
  long dims_data2[2] = {1, this->nslopes};
  this->d_centroids_ref = new carma_obj<T>(this->current_context, dims_data2);
  this->d_centroids_ref->reset();
}

template <class Tin, class T>
sutra_centroider_maskedPix<Tin, T>::~sutra_centroider_maskedPix() {}

template <class Tin, class T>
string sutra_centroider_maskedPix<Tin, T>::get_type() {
  return "maskedPix";
}

template <class Tin, class T>
int sutra_centroider_maskedPix<Tin, T>::get_cog(float *img, T *intensities,
                                                T *centroids, int nvalid,
                                                int npix, int ntot) {
  // TODO(Implement sutra_centroider_maskedPix<Tin, T>::get_cog)

  return get_maskedPix(img, intensities, centroids, this->d_validx->getData(),
                       this->d_validy->getData(), this->nvalid, ntot);
}

template <class Tin, class T>
int sutra_centroider_maskedPix<Tin, T>::get_maskedPix(T *img, T *intensities,
                                                      T *centroids,
                                                      int *subindx,
                                                      int *subindy, int nvalid,
                                                      int ns) {
  this->current_context->set_activeDevice(this->device, 1);

  fill_intensities(this->d_intensities->getData(), img, subindx, subindy, ns,
                   this->nslopes,
                   this->current_context->get_device(this->device));

  // T p_sum = reduce<T>(this->d_intensities->getData(), this->nslopes);
  this->d_intensities->reduceCub();

  getMaskedPix<T>(centroids, this->d_centroids_ref->getData(), img, subindx,
                  subindy, this->d_intensities->getOData(), ns, this->nslopes,
                  this->current_context->get_device(this->device));

  return EXIT_SUCCESS;
}

template <class Tin, class T>
int sutra_centroider_maskedPix<Tin, T>::get_cog(T *intensities, T *slopes,
                                                bool noise) {
  if (this->wfs != nullptr) {
    if (this->wfs->type == "pyrhr") {
      if (noise || this->wfs->roket == false) {
        return this->get_maskedPix(
            *(this->wfs->d_binimg), intensities, slopes,
            *(this->wfs->d_validsubsx), *(this->wfs->d_validsubsy),
            this->wfs->nvalid, this->wfs->nfft / this->wfs->nrebin);
      } else
        return this->get_maskedPix(
            *(this->wfs->d_binimg_notnoisy), intensities, slopes,
            *(this->wfs->d_validsubsx), *(this->wfs->d_validsubsy),
            this->wfs->nvalid, this->wfs->nfft / this->wfs->nrebin);
    } else
      DEBUG_TRACE("WFS must be pyrhr");
  }
  DEBUG_TRACE("this->wfs was not initialized");
  return EXIT_FAILURE;
}

template <class Tin, class T>
int sutra_centroider_maskedPix<Tin, T>::get_cog() {
  if (this->wfs != nullptr)
    return this->get_cog(*(this->wfs->d_intensities), *(this->wfs->d_slopes),
                         true);

  DEBUG_TRACE("this->wfs was not initialized");
  return EXIT_FAILURE;
}

template class sutra_centroider_maskedPix<float, float>;
template class sutra_centroider_maskedPix<uint16_t, float>;
