#include <sutra_centroider_maskedPix.h>
#include <string>

sutra_centroider_maskedPix::sutra_centroider_maskedPix(
    carma_context *context, sutra_wfs *wfs, long nvalid, long npupils,
    float offset, float scale, int device)
    : sutra_centroider(context, wfs, nvalid, offset, scale, device) {
  context->set_activeDevice(device, 1);

  if (wfs != nullptr)
    this->nslopes = wfs->d_validsubsx->getNbElem();
  else
    this->nslopes = nvalid * npupils;
  long dims_data[2] = {1, this->nslopes};
  this->d_intensities = new carma_obj<float>(context, dims_data);
  this->d_intensities->init_reduceCub();
  long dims_data2[2] = {1, nslopes};
  this->d_centroids_ref =
      new carma_obj<float>(this->current_context, dims_data2);
  this->d_centroids_ref->reset();
}

sutra_centroider_maskedPix::~sutra_centroider_maskedPix() {}

string sutra_centroider_maskedPix::get_type() { return "maskedPix"; }

int sutra_centroider_maskedPix::get_cog(float *img, float *intensities,
                                        float *centroids, int nvalid, int npix,
                                        int ntot) {
  // TODO(Implement sutra_centroider_maskedPix::get_cog)

  return get_maskedPix(img, intensities, centroids, this->d_validx->getData(),
                       this->d_validy->getData(), this->nvalid, ntot);
}

int sutra_centroider_maskedPix::get_maskedPix(float *img, float *intensities,
                                              float *centroids, int *subindx,
                                              int *subindy, int nvalid,
                                              int ns) {
  current_context->set_activeDevice(device, 1);

  fill_intensities(this->d_intensities->getData(), img, subindx, subindy, ns,
                   this->nslopes, this->current_context->get_device(device));

  // float p_sum = reduce<float>(this->d_intensities->getData(), this->nslopes);
  this->d_intensities->reduceCub();

  getMaskedPix<float>(centroids, this->d_centroids_ref->getData(), img, subindx,
                      subindy, this->d_intensities->getOData(), ns,
                      this->nslopes, this->current_context->get_device(device));

  return EXIT_SUCCESS;
}

int sutra_centroider_maskedPix::get_cog(float *intensities, float *slopes,
                                        bool noise) {
  if (this->wfs != nullptr) {
    if (wfs->type == "pyrhr") {
      if (noise || wfs->roket == false) {
        return this->get_maskedPix(*(wfs->d_binimg), intensities, slopes,
                                   *(wfs->d_validsubsx), *(wfs->d_validsubsy),
                                   wfs->nvalid, wfs->nfft / wfs->nrebin);
      } else
        return this->get_maskedPix(*(wfs->d_binimg_notnoisy), intensities,
                                   slopes, *(wfs->d_validsubsx),
                                   *(wfs->d_validsubsy), wfs->nvalid,
                                   wfs->nfft / wfs->nrebin);
    } else
      DEBUG_TRACE("WFS must be pyrhr");
  }
  DEBUG_TRACE("this->wfs was not initialized");
  return EXIT_FAILURE;
}

int sutra_centroider_maskedPix::get_cog() {
  if (this->wfs != nullptr)
    return this->get_cog(*(wfs->d_intensities), *(wfs->d_slopes), true);

  DEBUG_TRACE("this->wfs was not initialized");
  return EXIT_FAILURE;
}
