#include <sutra_centroider_bpcog.h>

sutra_centroider_bpcog::sutra_centroider_bpcog(carma_context *context,
    sutra_sensors *sensors, int nwfs, long nvalid, float offset, float scale, int device, int nmax) {

  this->current_context = context;

  this->device = device;
  context->set_activeDevice(device,1);
  this->wfs = sensors->d_wfs[nwfs];
  this->nwfs = nwfs;
  this->nvalid = nvalid;
  this->offset = offset;
  this->scale = scale;

  this->nmax = nmax;

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

string sutra_centroider_bpcog::get_type() {
  return "bpcog";
}

int sutra_centroider_bpcog::set_nmax(int nmax) {
  current_context->set_activeDevice(device,1);
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

int sutra_centroider_bpcog::get_cog(carma_streams *streams, float *cube,
                                    float *subsum, float *centroids, int nvalid, int npix, int ntot) {
  current_context->set_activeDevice(device,1);
  // brightest pixels cog
  // TODO: implemente sutra_centroider_bpcog::get_cog_async
  subap_sortmax<float>(npix * npix, nvalid, cube,this->d_bpix->getData(),this->d_bpind->getData(),this->nmax,current_context->get_device(device));
//carmaSafeCall(cudaDeviceSynchronize());
  /*
  	int nb = (int)(this->d_bpix->getNbElem());
  	float *tmpp;
  	tmpp=(float*)malloc((nb)*sizeof(float));
  	this->d_bpix->copyInto(tmpp,nb);
  	      for (int ii = 0 ; ii < nb ; ii++){
  	    	  printf("%5.5f \n",tmpp[ii]);
  	      }
  */
  subap_bpcentro<float>(this->nmax, nvalid, npix, this->d_bpix->getData(),this->d_bpind->getData(),
                        centroids, this->scale, this->offset);
  /*
  #if 1
    subap_centromax(npix * npix, nvalid, cube, centroids, npix, this->nmax,
        this->scale, this->offset);
  #else
    float *d_minim;
    cudaMalloc((void**)&d_minim, nvalid*sizeof(float));
    subap_centromax2<float>(npix * npix, nvalid, cube, centroids, d_minim, npix, this->nmax,
        this->scale, this->offset);
    cudaFree(d_minim);
  #endif
  */
  return EXIT_SUCCESS;
}

int sutra_centroider_bpcog::get_cog(float *subsum, float *slopes, bool noise) {
  if(noise || wfs->error_budget == false) {
    return this->get_cog(wfs->streams, *wfs->d_bincube, subsum, slopes,
                         wfs->nvalid, wfs->npix, wfs->d_bincube->getNbElem());
  } else {
    return this->get_cog(wfs->streams, *wfs->d_bincube_notnoisy, subsum, slopes,
                         wfs->nvalid, wfs->npix, wfs->d_bincube->getNbElem());
  }
}

int sutra_centroider_bpcog::get_cog() {
  return this->get_cog(*(wfs->d_subsum),*(wfs->d_slopes),true);
}
