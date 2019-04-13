#include <sutra_centroider_pyr.h>
#include <iostream>
#include <string>

template <class Tin, class T>
sutra_centroider_pyr<Tin, T>::sutra_centroider_pyr(carma_context *context,
                                                   sutra_wfs *wfs, long nvalid,
                                                   float offset, float scale,
                                                   bool filter_TT, int device)
  : sutra_centroider<Tin, T>(context, wfs, nvalid, offset, scale, filter_TT, device) {
  context->set_activeDevice(device, 1);

  this->nslopes = 2 * nvalid;

  this->pyr_type = "pyrhr";

  this->valid_thresh = 1e-4;

  // centroider method by default nosin_global
  this->method = Method_CoG(false, false);

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
sutra_centroider_pyr<Tin, T>::~sutra_centroider_pyr() {}

template <class Tin, class T>
string sutra_centroider_pyr<Tin, T>::get_type() {
  return this->pyr_type;
}

template <class Tin, class T>
int sutra_centroider_pyr<Tin, T>::set_valid_thresh(T valid_thresh) {
  this->valid_thresh = valid_thresh;
  return EXIT_SUCCESS;
}
template <class Tin, class T>
T sutra_centroider_pyr<Tin, T>::get_valid_thresh() {
  return this->valid_thresh;
}

template <class Tin, class T>
int sutra_centroider_pyr<Tin, T>::set_method(Method_CoG method) {
  this->method = method;
  return EXIT_SUCCESS;
}

template <class Tin, class T>
Method_CoG sutra_centroider_pyr<Tin, T>::get_method() {
  return this->method;
}

template <class Tin, class T>
string sutra_centroider_pyr<Tin, T>::get_method_str() {
  return Method_CoG::str(this->method);
}

template <class Tin, class T>
int sutra_centroider_pyr<Tin, T>::get_cog(float *cube, float *intensities,
                                          T *centroids, int nvalid, int npix,
                                          int ntot) {
  // TODO(Implement sutra_centroider_pyr<Tin, T>::get_cog)

  return get_pyr(cube, intensities, centroids, this->d_validx->getData(),
                 this->d_validy->getData(), this->nvalid,
                 this->d_bincube->getDims(1), 4);
}

template <class Tin, class T>
int sutra_centroider_pyr<Tin, T>::get_pyr(float *cube, float *intensities,
                                          T *centroids, int *subindx,
                                          int *subindy, int nvalid, int ns,
                                          int nim) {
  this->current_context->set_activeDevice(this->device, 1);

  pyr_intensities(intensities, cube, subindx, subindy, ns, nvalid, nim,
                  this->current_context->get_device(this->device));

  if (!(this->method.isLocal)) {
    float p_sum = reduce<float>(intensities, nvalid);
    fillvalues<float>(intensities, p_sum, nvalid,
               this->current_context->get_device(this->device));
  }

  pyr2_slopes(centroids, this->d_centroids_ref->getData(), cube, subindx,
              subindy, intensities, ns, nvalid, this->scale, this->valid_thresh,
              this->method.isSinus,  // if we are using a sin method
              this->current_context->get_device(this->device));

  carma_axpy<float>(this->current_context->get_cublasHandle(), this->nslopes,
                    -1.0f, this->d_centroids_ref->getData(), 1, centroids, 1);

  return EXIT_SUCCESS;
}

template <class Tin, class T>
int sutra_centroider_pyr<Tin, T>::get_cog(float *intensities, T *slopes,
                                          bool noise) {
  if (this->wfs != nullptr) {
    if (this->pyr_type == "pyr" || this->pyr_type == "roof")
      return this->get_pyr(*(this->wfs->d_bincube), intensities, slopes,
                           *(this->wfs->d_validsubsx),
                           *(this->wfs->d_validsubsy), this->wfs->nvalid,
                           this->wfs->nfft / this->wfs->nrebin, 4);
    else if (this->pyr_type == "pyrhr") {
      if (noise || this->wfs->roket == false) {
        return this->get_pyr(*(this->wfs->d_binimg), intensities, slopes,
                             *(this->wfs->d_validsubsx),
                             *(this->wfs->d_validsubsy), this->wfs->nvalid,
                             this->wfs->nfft / this->wfs->nrebin, 4);
      } else
        return this->get_pyr(*(this->wfs->d_binimg_notnoisy), intensities,
                             slopes, *(this->wfs->d_validsubsx),
                             *(this->wfs->d_validsubsy), this->wfs->nvalid,
                             this->wfs->nfft / this->wfs->nrebin, 4);
    }
  }

  DEBUG_TRACE("this->wfs was not initialized");
  return EXIT_FAILURE;
}

template <class Tin, class T>
int sutra_centroider_pyr<Tin, T>::get_cog() {
  if (this->wfs != nullptr)
    return this->get_cog(*(this->wfs->d_intensities), *(this->wfs->d_slopes),
                         true);

  DEBUG_TRACE("this->wfs was not initialized");
  return EXIT_FAILURE;
}

template class sutra_centroider_pyr<float, float>;
template class sutra_centroider_pyr<uint16_t, float>;
