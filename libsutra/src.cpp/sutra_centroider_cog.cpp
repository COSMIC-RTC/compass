#include <sutra_centroider_cog.h>
#include <string>

sutra_centroider_cog::sutra_centroider_cog(carma_context *context, sutra_sensors *sensors, int nwfs,
    long nvalid, float offset, float scale, int device) {
  this->current_context = context;

  this->device = device;
  context->set_activeDevice(device,1);
  this->wfs = sensors->d_wfs[nwfs];
  this->nwfs = nwfs;
  this->nvalid = nvalid;
  this->offset = offset;
  this->scale = scale;

}

sutra_centroider_cog::~sutra_centroider_cog() {
}

string sutra_centroider_cog::get_type() {
  return "cog";
}

int sutra_centroider_cog::get_cog(carma_streams *streams, float *cube,
    float *subsum, float *centroids, int nvalid, int npix, int ntot) {

  current_context->set_activeDevice(device,1);
  // simple cog
  int nstreams = streams->get_nbStreams();
  //fprintf(stderr, "\n[%s@%d]: nstreams=%d\n", __FILE__, __LINE__, nstreams);
  if (nstreams > 1) {
    //fprintf(stderr, "\n[%s@%d]: i=%d istart=%d npix=%d nvalid=%d\n", __FILE__, __LINE__, i, istart, npix, nvalid);
    subap_reduce_async(npix * npix, nvalid, streams, cube, subsum);
    //fprintf(stderr, "\n[%s@%d] I'm here\n", __FILE__, __LINE__);
    get_centroids_async(npix * npix, nvalid, npix, streams, cube, centroids,
        subsum, this->scale, this->offset);
    //fprintf(stderr, "\n[%s@%d] I'm here\n", __FILE__, __LINE__);
    //streams->wait_all_streams();
  } else {
    subap_reduce(ntot, (npix * npix), nvalid, cube, subsum, current_context->get_device(device));
    get_centroids(ntot, (npix * npix), nvalid, npix, cube, centroids, subsum,
        this->scale, this->offset, current_context->get_device(device));
  }
  return EXIT_SUCCESS;
}

int sutra_centroider_cog::get_cog(float *subsum, float *slopes) {
  return this->get_cog(wfs->streams, *(wfs->d_bincube), subsum,
      slopes, wfs->nvalid_tot, wfs->npix, wfs->d_bincube->getNbElem());
}

int sutra_centroider_cog::get_cog() {
  return this->get_cog(*(wfs->d_subsum),*(wfs->d_slopes));
}

