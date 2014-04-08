#include <sutra_centroider_pyr.h>
#include <string>

sutra_centroider_pyr::sutra_centroider_pyr(carma_context *context, long nwfs,
    long nvalid, float offset, float scale, int device) {
  this->current_context = context;

  this->device = device;
  context->set_activeDevice(device);
  this->nwfs = nwfs;
  this->nvalid = nvalid;
  this->offset = offset;
  this->scale = scale;

}

sutra_centroider_pyr::~sutra_centroider_pyr() {

}

string sutra_centroider_pyr::get_type() {
  return "pyr";
}

int sutra_centroider_pyr::init_bincube(sutra_wfs *wfs) {

  return EXIT_SUCCESS;
}

int sutra_centroider_pyr::get_cog(carma_streams *streams, float *cube,
    float *subsum, float *centroids, int nvalid, int npix, int ntot) {
  //TODO: Implement sutra_centroider_pyr::get_cog
  cerr << "get_cog not implemented\n";

  return EXIT_SUCCESS;
}

int sutra_centroider_pyr::get_pyr(float *cube, float *subsum, float *centroids,
    int *subindx, int *subindy, int nvalid, int ns, int nim) {

  pyr_slopes(centroids, cube, subindx, subindy, subsum, ns, nvalid, nim,
      this->device);

  return EXIT_SUCCESS;
}

int sutra_centroider_pyr::get_cog(sutra_wfs *wfs, float *slopes) {
  return this->get_pyr(*(wfs->d_bincube), *(wfs->d_subsum), slopes,
      *(wfs->d_validsubsx), *(wfs->d_validsubsy), wfs->nvalid,
      wfs->nfft / wfs->nrebin, 4);
}

int sutra_centroider_pyr::get_cog(sutra_wfs *wfs) {
  return this->get_cog(wfs, *(wfs->d_slopes));
}

