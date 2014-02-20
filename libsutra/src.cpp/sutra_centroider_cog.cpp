#include <sutra_centroider_cog.h>
#include <string>

sutra_centroider_cog::sutra_centroider_cog(carma_context *context, long nwfs,
    long nvalid, float offset, float scale, int device) {
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

string sutra_centroider_cog::get_type(){
  return "cog";
}

int sutra_centroider_cog::init_bincube(sutra_wfs *wfs) {

  return EXIT_SUCCESS;
}

int sutra_centroider_cog::get_cog(carma_streams *streams, float *cube,
    float *subsum, float *centroids, int nvalid, int npix, int ntot) {
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
    subap_reduce(ntot, npix * npix, nvalid, cube, subsum);

    get_centroids(ntot, npix * npix, nvalid, npix, cube, centroids, subsum,
        this->scale, this->offset);
  }
  return EXIT_SUCCESS;
}

int sutra_centroider_cog::get_cog(sutra_wfs *wfs, float *slopes) {
  return this->get_cog(wfs->streams, *(wfs->d_bincube), *(wfs->d_subsum),
      slopes, wfs->nvalid, wfs->npix, wfs->d_bincube->getNbElem());
}

int sutra_centroider_cog::get_cog(sutra_wfs *wfs) {
  return this->get_cog(wfs, *(wfs->d_slopes));
}

