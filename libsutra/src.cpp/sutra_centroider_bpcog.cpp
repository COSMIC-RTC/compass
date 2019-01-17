#include <sutra_centroider_bpcog.h>

sutra_centroider_bpcog::sutra_centroider_bpcog(carma_context *context,
                                               sutra_wfs *wfs, long nvalid,
                                               float offset, float scale,
                                               int device, int nmax)
    : sutra_centroider(context, wfs, nvalid, offset, scale, device) {
  this->nslopes = 2 * nvalid;
  this->nmax = nmax;
  long dims_data2[2] = {1, nslopes};
  this->d_centroids_ref =
      new carma_obj<float>(this->current_context, dims_data2);
  this->d_centroids_ref->reset();

  long dims_data[3];
  dims_data[0] = 2;
  dims_data[1] = nmax;
  dims_data[2] = nvalid;

  this->d_bpix = new carma_obj<float>(this->current_context, dims_data);
  this->d_bpind = new carma_obj<uint>(this->current_context, dims_data);
}

sutra_centroider_bpcog::~sutra_centroider_bpcog() {
  delete this->d_bpix;
  delete this->d_bpind;
}

string sutra_centroider_bpcog::get_type() { return "bpcog"; }

int sutra_centroider_bpcog::set_nmax(int nmax) {
  current_context->set_activeDevice(device, 1);
  this->nmax = nmax;
  delete this->d_bpix;
  delete this->d_bpind;

  long dims_data[3];
  dims_data[0] = 2;
  dims_data[1] = nmax;
  dims_data[2] = this->nvalid;

  this->d_bpix = new carma_obj<float>(this->current_context, dims_data);
  this->d_bpind = new carma_obj<uint>(this->current_context, dims_data);

  return EXIT_SUCCESS;
}

int sutra_centroider_bpcog::get_cog(float *img, float *intensities,
                                    float *centroids, int nvalid, int npix,
                                    int ntot) {
  current_context->set_activeDevice(device, 1);

  get_centroids(ntot, (npix * npix), nvalid, npix, img, centroids,
                this->d_centroids_ref->getData(), this->d_validx->getData(),
                this->d_validy->getData(), intensities, this->nmax, this->scale,
                this->offset, current_context->get_device(device));

  // brightest pixels cog
  // subap_sortmax<float>(npix * npix, nvalid, cube, this->d_bpix->getData(),
  //                      this->d_bpind->getData(), this->nmax,
  //                      current_context->get_device(device));

  // subap_bpcentro<float>(this->nmax, nvalid, npix, this->d_bpix->getData(),
  //                       this->d_bpind->getData(), centroids, this->scale,
  //                       this->offset);

  return EXIT_SUCCESS;
}

int sutra_centroider_bpcog::get_cog(float *intensities, float *slopes,
                                    bool noise) {
  if (this->wfs != nullptr) {
    if (noise || wfs->roket == false)
      return this->get_cog(*(wfs->d_binimg), intensities, slopes, wfs->nvalid,
                           wfs->npix, wfs->d_binimg->getDims()[1]);
    else
      return this->get_cog(*(wfs->d_binimg_notnoisy), intensities, slopes,
                           wfs->nvalid, wfs->npix, wfs->d_binimg->getDims()[1]);
  }
  DEBUG_TRACE("this->wfs was not initialized");
  return EXIT_FAILURE;
}

int sutra_centroider_bpcog::get_cog() {
  if (this->wfs != nullptr)
    return this->get_cog(*(wfs->d_intensities), *(wfs->d_slopes), true);

  DEBUG_TRACE("this->wfs was not initialized");
  return EXIT_FAILURE;
}
