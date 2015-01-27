#include <sutra_centroider_tcog.h>
#include <string>

sutra_centroider_tcog::sutra_centroider_tcog(carma_context *context, long nwfs,
    long nvalid, float offset, float scale, int device) {
  this->current_context = context;

  this->nwfs = nwfs;
  this->nvalid = nvalid;
  this->device = device;
  context->set_activeDevice(device);
  this->offset = offset;
  this->scale = scale;

  this->threshold = 0;

}

sutra_centroider_tcog::~sutra_centroider_tcog() {

}

string sutra_centroider_tcog::get_type() {
  return "tcog";
}

int sutra_centroider_tcog::init_bincube(sutra_wfs *wfs) {
  return EXIT_SUCCESS;
}

int sutra_centroider_tcog::set_threshold(float threshold) {
  this->threshold = threshold;

  return EXIT_SUCCESS;
}

int sutra_centroider_tcog::get_cog(carma_streams *streams, float *cube,
    float *subsum, float *centroids, int nvalid, int npix, int ntot) {
  //TODO: Implement sutra_centroider_tcog::get_cog_async
  subap_reduce(ntot, npix * npix, nvalid, cube, subsum, this->threshold, this->current_context->get_device(device));

  get_centroids(ntot, npix * npix, nvalid, npix, cube, centroids, subsum,
      this->threshold, this->scale, this->offset, this->current_context->get_device(device));

  return EXIT_SUCCESS;
}

int sutra_centroider_tcog::get_cog(sutra_wfs *wfs, float *slopes) {
  return this->get_cog(wfs->streams, *(wfs->d_bincube), *(wfs->d_subsum),
      slopes, wfs->nvalid, wfs->npix, wfs->d_bincube->getNbElem());
}

int sutra_centroider_tcog::get_cog(sutra_wfs *wfs) {
  return this->get_cog(wfs, *(wfs->d_slopes));
}

