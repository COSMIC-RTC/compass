#include <sutra_centroider_maskedPix.h>
#include <string>

sutra_centroider_maskedPix::sutra_centroider_maskedPix(
    carma_context *context, sutra_wfs *wfs, long nvalid, long npupils,
    float offset, float scale, int device)
    : sutra_centroider(context, wfs, nvalid, offset, scale, device) {
  context->set_activeDevice(device, 1);

  this->nslopes = wfs->d_validsubsx->getNbElem();
  long dims_data[2] = {1, this->nslopes};
  this->d_subsum = new carma_obj<float>(context, dims_data);
}

sutra_centroider_maskedPix::~sutra_centroider_maskedPix() {}

string sutra_centroider_maskedPix::get_type() { return "maskedPix"; }

int sutra_centroider_maskedPix::get_cog(carma_streams *streams, float *cube,
                                        float *subsum, float *centroids,
                                        int nvalid, int npix, int ntot) {
  // TODO(Implement sutra_centroider_maskedPix::get_cog)

  return get_maskedPix(cube, subsum, centroids, this->d_validx->getData(),
                       this->d_validy->getData(), this->nvalid,
                       this->d_bincube->getDims(1), 4);
}

int sutra_centroider_maskedPix::get_maskedPix(float *cube, float *subsum,
                                              float *centroids, int *subindx,
                                              int *subindy, int nvalid, int ns,
                                              int nim) {
  current_context->set_activeDevice(device, 1);

  fill_subsum(this->d_subsum->getData(), cube, subindx, subindy, ns,
              this->nslopes, this->current_context->get_device(device));

  float p_sum = reduce<float>(this->d_subsum->getData(), this->nslopes);

  getMaskedPix<float>(centroids, cube, subindx, subindy, p_sum, ns,
                      this->nslopes, this->current_context->get_device(device));

  return EXIT_SUCCESS;
}

int sutra_centroider_maskedPix::get_cog(float *subsum, float *slopes,
                                        bool noise) {
  if (this->wfs != nullptr) {
    if (wfs->type == "pyrhr") {
      if (noise || wfs->roket == false) {
        return this->get_maskedPix(*(wfs->d_binimg), subsum, slopes,
                                   *(wfs->d_validsubsx), *(wfs->d_validsubsy),
                                   wfs->nvalid, wfs->nfft / wfs->nrebin, 4);
      } else
        return this->get_maskedPix(*(wfs->d_binimg_notnoisy), subsum, slopes,
                                   *(wfs->d_validsubsx), *(wfs->d_validsubsy),
                                   wfs->nvalid, wfs->nfft / wfs->nrebin, 4);
    } else
      DEBUG_TRACE("WFS must be pyrhr");
  }
  DEBUG_TRACE("this->wfs was not initialized");
  return EXIT_FAILURE;
}

int sutra_centroider_maskedPix::get_cog() {
  if (this->wfs != nullptr)
    return this->get_cog(*(wfs->d_subsum), *(wfs->d_slopes), true);

  DEBUG_TRACE("this->wfs was not initialized");
  return EXIT_FAILURE;
}
