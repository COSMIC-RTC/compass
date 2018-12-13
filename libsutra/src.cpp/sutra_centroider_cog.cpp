#include <sutra_centroider_cog.h>
#include <string>

sutra_centroider_cog::sutra_centroider_cog(carma_context *context,
                                           sutra_wfs *wfs, long nvalid,
                                           float offset, float scale,
                                           int device)
    : sutra_centroider(context, wfs, nvalid, offset, scale, device) {
  this->nslopes = 2 * nvalid;
  long dims_data2[2] = {1, nslopes};
  this->d_centroids_ref =
      new carma_obj<float>(this->current_context, dims_data2);
}

sutra_centroider_cog::~sutra_centroider_cog() {}

string sutra_centroider_cog::get_type() { return "cog"; }

int sutra_centroider_cog::get_cog(float *img, float *intensities,
                                  float *centroids, int nvalid, int npix,
                                  int ntot) {
  current_context->set_activeDevice(device, 1);

  // subap_reduce(ntot, (npix * npix), nvalid, img, ref,
  //              current_context->get_device(device));
  get_centroids(ntot, (npix * npix), nvalid, npix, img, centroids,
                this->d_centroids_ref->getData(), this->d_validx->getData(),
                this->d_validy->getData(), intensities, this->scale,
                this->offset, current_context->get_device(device));

  return EXIT_SUCCESS;
}

int sutra_centroider_cog::get_cog(float *intensities, float *slopes,
                                  bool noise) {
  if (this->wfs != nullptr) {
    if (noise || wfs->roket == false) {
      return this->get_cog(*(wfs->d_binimg), intensities, slopes,
                           wfs->nvalid_tot, wfs->npix,
                           wfs->d_binimg->getDims()[1]);
    } else {
      return this->get_cog(*(wfs->d_binimg_notnoisy), intensities, slopes,
                           wfs->nvalid_tot, wfs->npix,
                           wfs->d_binimg->getDims()[1]);
    }
  }
  DEBUG_TRACE("this->wfs was not initialized");
  return EXIT_FAILURE;
}

int sutra_centroider_cog::get_cog() {
  if (this->wfs != nullptr)
    return this->get_cog(*(wfs->d_intensities), *(wfs->d_slopes), true);

  DEBUG_TRACE("this->wfs was not initialized");
  return EXIT_FAILURE;
}
