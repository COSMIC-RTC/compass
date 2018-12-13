#include <sutra_centroider_pyr.h>
#include <iostream>
#include <string>

sutra_centroider_pyr::sutra_centroider_pyr(carma_context *context,
                                           sutra_wfs *wfs, long nvalid,
                                           float offset, float scale,
                                           int device)
    : sutra_centroider(context, wfs, nvalid, offset, scale, device) {
  context->set_activeDevice(device, 1);

  this->nslopes = 2 * nvalid;

  this->pyr_type = "pyrhr";

  this->valid_thresh = 1e-4;

  // centroider method by default nosin_global
  this->method = Method_CoG(false, false);
}

sutra_centroider_pyr::~sutra_centroider_pyr() {}

string sutra_centroider_pyr::get_type() { return this->pyr_type; }

int sutra_centroider_pyr::set_valid_thresh(float valid_thresh) {
  this->valid_thresh = valid_thresh;
  return EXIT_SUCCESS;
}
float sutra_centroider_pyr::get_valid_thresh() { return this->valid_thresh; }

int sutra_centroider_pyr::set_method(Method_CoG method) {
  this->method = method;
  return EXIT_SUCCESS;
}

Method_CoG sutra_centroider_pyr::get_method() { return this->method; }

string sutra_centroider_pyr::get_method_str() {
  return Method_CoG::str(this->method);
}

int sutra_centroider_pyr::get_cog(float *cube, float *intensities,
                                  float *centroids, int nvalid, int npix,
                                  int ntot) {
  // TODO(Implement sutra_centroider_pyr::get_cog)

  return get_pyr(cube, intensities, centroids, this->d_validx->getData(),
                 this->d_validy->getData(), this->nvalid,
                 this->d_bincube->getDims(1), 4);
}

int sutra_centroider_pyr::get_pyr(float *cube, float *intensities,
                                  float *centroids, int *subindx, int *subindy,
                                  int nvalid, int ns, int nim) {
  current_context->set_activeDevice(device, 1);

  pyr_intensities(intensities, cube, subindx, subindy, ns, nvalid, nim,
                  this->current_context->get_device(device));

  if (!(this->method.isLocal)) {
    float p_sum = reduce<float>(intensities, nvalid);
    fillvalues(intensities, p_sum, nvalid,
               this->current_context->get_device(device));
  }

  pyr2_slopes(centroids, cube, subindx, subindy, intensities, ns, nvalid,
              this->scale, this->valid_thresh,
              this->method.isSinus,  // if we are using a sin method
              this->current_context->get_device(device));

  return EXIT_SUCCESS;
}

int sutra_centroider_pyr::get_cog(float *intensities, float *slopes,
                                  bool noise) {
  if (this->wfs != nullptr) {
    if (this->pyr_type == "pyr" || this->pyr_type == "roof")
      return this->get_pyr(*(wfs->d_bincube), intensities, slopes,
                           *(wfs->d_validsubsx), *(wfs->d_validsubsy),
                           wfs->nvalid, wfs->nfft / wfs->nrebin, 4);
    else if (this->pyr_type == "pyrhr") {
      if (noise || wfs->roket == false) {
        return this->get_pyr(*(wfs->d_binimg), intensities, slopes,
                             *(wfs->d_validsubsx), *(wfs->d_validsubsy),
                             wfs->nvalid, wfs->nfft / wfs->nrebin, 4);
      } else
        return this->get_pyr(*(wfs->d_binimg_notnoisy), intensities, slopes,
                             *(wfs->d_validsubsx), *(wfs->d_validsubsy),
                             wfs->nvalid, wfs->nfft / wfs->nrebin, 4);
    }
  }

  DEBUG_TRACE("this->wfs was not initialized");
  return EXIT_FAILURE;
}

int sutra_centroider_pyr::get_cog() {
  if (this->wfs != nullptr)
    return this->get_cog(*(wfs->d_intensities), *(wfs->d_slopes), true);

  DEBUG_TRACE("this->wfs was not initialized");
  return EXIT_FAILURE;
}
