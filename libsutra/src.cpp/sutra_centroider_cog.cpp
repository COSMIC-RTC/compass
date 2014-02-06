#include <sutra_centroider_cog.h>
#include <string>

sutra_centroider_cog::sutra_centroider_cog(carma_context *context, long nwfs,
    long nvalid, float offset, float scale, int device) {
  this->type = "cog";

  this->current_context = context;

  this->device = device;
  context->set_activeDevice(device);
  this->nwfs = nwfs;
  this->nvalid = nvalid;
  this->offset = offset;
  this->scale = scale;

}

sutra_centroider_cog::~sutra_centroider_cog() {
}

int sutra_centroider_cog::init_bincube(sutra_wfs *wfs) {

  return EXIT_SUCCESS;
}

int sutra_centroider_cog::get_cog(float *cube, float *subsum, float *centroids,
    int nvalid, int npix, int ntot) {
  // simple cog

  subap_reduce(ntot, npix * npix, nvalid, cube, subsum);

  get_centroids(ntot, npix * npix, nvalid, npix, cube, centroids, subsum,
      this->scale, this->offset);

  return EXIT_SUCCESS;
}

int sutra_centroider_cog::get_cog(sutra_wfs *wfs, carma_obj<float> *slopes) {
  return this->get_cog(wfs->d_bincube->getData(), wfs->d_subsum->getData(),
      slopes->getData(), wfs->nvalid, wfs->npix, wfs->d_bincube->getNbElem());
}

int sutra_centroider_cog::get_cog(sutra_wfs *wfs) {
  return this->get_cog(wfs, wfs->d_slopes);
}

int sutra_centroider_cog::get_cog_async(carma_streams *streams, float *cube,
    float *subsum, float *centroids, int nvalid, int npix) {
  // simple cog
  subap_reduce_async(npix * npix, nvalid, streams, cube, subsum);
  get_centroids_async(npix * npix, nvalid, npix, streams, cube, centroids,
      subsum, this->scale, this->offset);
  return EXIT_SUCCESS;
}

int sutra_centroider_cog::get_cog_async(sutra_wfs *wfs,
    carma_obj<float> *slopes) {
  return this->get_cog_async(wfs->streams, wfs->d_bincube->getData(),
      wfs->d_subsum->getData(), slopes->getData(), wfs->nvalid, wfs->npix);
}

int sutra_centroider_cog::get_cog_async(sutra_wfs *wfs) {
  return this->get_cog_async(wfs, wfs->d_slopes);
}

