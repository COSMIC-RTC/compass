#include <sutra_centroider_bpcog.h>

sutra_centroider_bpcog::sutra_centroider_bpcog(carma_context *context,
    long nwfs, long nvalid, float offset, float scale, int device, int nmax) {

  this->current_context = context;

  this->device = device;
  context->set_activeDevice(device);
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

int sutra_centroider_bpcog::init_bincube(sutra_wfs *wfs) {
  return EXIT_SUCCESS;
}

int sutra_centroider_bpcog::set_nmax(int nmax) {
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
  // brightest pixels cog
  // TODO: implemente sutra_centroider_bpcog::get_cog_async
	subap_sortmax<float>(npix * npix, nvalid, cube,this->d_bpix->getData(),this->d_bpind->getData(),this->nmax,current_context->get_device(device));
//cutilSafeCall(cudaDeviceSynchronize());
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

int sutra_centroider_bpcog::get_cog(sutra_wfs *wfs, float *slopes) {
  return this->get_cog(wfs->streams, *wfs->d_bincube, *wfs->d_subsum, slopes,
      wfs->nvalid, wfs->npix, wfs->d_bincube->getNbElem());
}

int sutra_centroider_bpcog::get_cog(sutra_wfs *wfs) {
  return this->get_cog(wfs, *wfs->d_slopes);
}
