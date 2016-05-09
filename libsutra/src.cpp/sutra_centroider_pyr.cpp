#include <sutra_centroider_pyr.h>
#include <sutra_wfs_pyr_pyr4.h>
#include <string>

sutra_centroider_pyr::sutra_centroider_pyr(carma_context *context,
                                           sutra_sensors *sensors, int nwfs,
                                           long nvalid, float offset,
                                           float scale, int device) {
  if (sensors->d_wfs[nwfs]->type != "pyr")
    throw "sutra_centroider_roof expect a sutra_wfs_pyr_pyr4 sensor";
  this->current_context = context;

  this->device = device;
  context->set_activeDevice(device,1);
  this->wfs = sensors->d_wfs[nwfs];
  this->nwfs = nwfs;
  this->nvalid = nvalid;
  this->offset = offset;
  this->scale = scale;
  //this->pyr_type = pyr_type;

}

sutra_centroider_pyr::~sutra_centroider_pyr() {

}

string sutra_centroider_pyr::get_type() {
  return "pyr";
}

int sutra_centroider_pyr::get_cog(carma_streams *streams, float *cube,
                                  float *subsum, float *centroids, int nvalid,
                                  int npix, int ntot) {
  //TODO: Implement sutra_centroider_pyr::get_cog
  cerr << "get_cog not implemented\n";

  return EXIT_SUCCESS;
}

int sutra_centroider_pyr::get_pyr(float *cube, float *subsum, float *centroids,
                                  int *subindx, int *subindy, int nvalid,
                                  int ns, int nim) {
  current_context->set_activeDevice(device,1);
  //pyr_slopes(centroids, cube, subindx, subindy, subsum, ns, nvalid, nim, this->current_context->get_device(device));
  pyr2_slopes(centroids, cube, subindx, subindy, subsum, ns, nvalid,
              this->current_context->get_device(device));
  return EXIT_SUCCESS;
}

int sutra_centroider_pyr::get_cog(float *subsum, float *slopes) {
  /*
   return this->get_pyr(*(wfs->d_bincube), *(wfs->d_subsum), slopes,
   *(wfs->d_validsubsx), *(wfs->d_validsubsy), wfs->nvalid,
   wfs->nfft / wfs->nrebin, 4);
   */
  return this->get_pyr(*(wfs->d_binimg), *(wfs->d_subsum), slopes,
                       *(wfs->d_validsubsx), *(wfs->d_validsubsy), wfs->nvalid,
                       wfs->nfft / wfs->nrebin, 4);
}

int sutra_centroider_pyr::get_cog() {
  return this->get_cog(*(wfs->d_subsum), *(wfs->d_slopes));
}

