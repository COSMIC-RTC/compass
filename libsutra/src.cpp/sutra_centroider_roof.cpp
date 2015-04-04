#include <sutra_centroider_roof.h>
#include <string>

sutra_centroider_roof::sutra_centroider_roof(carma_context *context, sutra_sensors *sensors, int nwfs,
    long nvalid, float offset, float scale, int device) {
  if(sensors->d_wfs[nwfs]->type != "roof")
    throw "sutra_centroider_roof expect a sutra_wfs_pyr_roof sensor";
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

sutra_centroider_roof::~sutra_centroider_roof() {

}

string sutra_centroider_roof::get_type() {
  return "roof";
}

int sutra_centroider_roof::get_cog(carma_streams *streams, float *cube,
    float *subsum, float *centroids, int nvalid, int npix, int ntot) {
  //TODO: Implement sutra_centroider_roof::get_cog
  cerr << "get_cog not implemented\n";

  return EXIT_SUCCESS;
}

int sutra_centroider_roof::get_roof(float *cube, float *subsum,
    float *centroids, int *subindx, int *subindy, int nvalid, int ns, int nim) {
  current_context->set_activeDevice(device,1);
  roof_slopes(centroids, cube, subindx, subindy, subsum, ns, nvalid, nim,
      this->current_context->get_device(device));
  return EXIT_SUCCESS;
}

int sutra_centroider_roof::get_cog(float *subsum, float *slopes) {
  return this->get_roof(*(wfs->d_bincube), subsum, slopes,
    *(wfs->d_validsubsx), *(wfs->d_validsubsy), wfs->nvalid,
    wfs->nfft / wfs->nrebin, 4);
}

int sutra_centroider_roof::get_cog() {
  return this->get_cog(*(wfs->d_subsum),*(wfs->d_slopes));
}

