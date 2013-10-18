#include <yoga_centroider.h>
#include <string>


yoga_centroider::yoga_centroider(yoga_context *context, long nwfs, long nvalid, float offset, float scale, int device,const char *typec)
{
  this->d_bincube  = 0L;
  this->d_subsum  = 0L;
  this->d_refcentro = 0L;
  this->d_validpix = 0L;
  this->d_validindx = 0L;

  this->d_centroids = 0L;
  this->d_weights = 0L;

  this->d_corrfnct  = 0L;
  this->d_corrspot  = 0L;
  this->d_corrnorm  = 0L;
  this->d_corrmax   = 0L;
  this->d_corr      = 0L;
  this->d_interpmat = 0L;

  this->current_context=context;

  long *dims_data1 = new long[2]; dims_data1[0] = 1;
  //long *dims_data2 = new long[3]; dims_data2[0] = 2; 
  //long *dims_data3 = new long[4]; dims_data3[0] = 3;

  this->nwfs        = nwfs;
  this->nvalid      = nvalid;
  this->device      = device; context->set_activeDevice(device);
  this->typec       = string(typec);
  this->offset      = offset;
  this->scale       = scale;

  dims_data1[1]     = 2*nvalid;
  //this->d_centroids = new yoga_obj<float>(context,dims_data1);

  delete[] dims_data1;
}


yoga_centroider::~yoga_centroider()
{
  if (this->d_centroids  != 0L) delete this->d_centroids;
  //delete this->current_context;
  if (this->d_bincube  != 0L) delete this->d_bincube;
  if (this->d_subsum  != 0L) delete this->d_subsum;
  if (this->d_centroids != 0L) delete this->d_centroids;
  if (this->d_refcentro != 0L) delete this->d_refcentro;
  if (this->d_validpix != 0L) delete this->d_validpix;
  if (this->d_validindx != 0L) delete this->d_validindx;

  if (this->d_weights != 0L) delete this->d_weights;

  if (this->d_corrfnct  != 0L) delete this->d_corrfnct;
  if (this->d_corrspot  != 0L) delete this->d_corrspot;
  if (this->d_corrnorm  != 0L) delete this->d_corrnorm;
  if (this->d_corrmax   != 0L) delete this->d_corrmax;
  if (this->d_corr      != 0L) delete this->d_corr;
  if (this->d_interpmat != 0L) delete this->d_interpmat;
}

int yoga_centroider::init_bincube(yoga_wfs *wfs)
{
  if (this->d_bincube != 0L) delete this->d_bincube;

  this->npix = wfs->npix;

  long *dims_data3 = new long[4];
  dims_data3[0] = 3; 
  dims_data3[1] = this->npix; dims_data3[2] = this->npix; dims_data3[3] = this->nvalid; 

  current_context->set_activeDevice(device);
  this->d_bincube = new yoga_obj<float>(current_context,dims_data3);

  if (this->d_subsum != 0L) delete this->d_subsum;

  long *dims_data1 = new long[2]; dims_data1[0] = 1;
  dims_data1[1]    = this->nvalid;
  this->d_subsum   = new yoga_obj<float>(current_context,dims_data1);

  delete[] dims_data1;
  delete[] dims_data3;

 return EXIT_SUCCESS;
}

int yoga_centroider::init_weights(yoga_wfs *wfs)
{
  if (this->d_weights != 0L) delete this->d_weights;

  this->npix = wfs->npix;

  long *dims_data3 = new long[4];
  dims_data3[0] = 3; 
  dims_data3[1] = this->npix; dims_data3[2] = this->npix; dims_data3[3] = this->nvalid; 

  current_context->set_activeDevice(device);
  this->d_weights = new yoga_obj<float>(current_context,dims_data3);

  delete[] dims_data3;

  return EXIT_SUCCESS;
}

int yoga_centroider::set_threshold(float threshold)
{
  this->threshold = threshold;

  return EXIT_SUCCESS;
}

int yoga_centroider::set_nmax(int nmax)
{
  this->nmax = nmax;

  return EXIT_SUCCESS;
}

int yoga_centroider::load_weights(float *weights, int ndim)
{
  if (ndim == 3) this->d_weights->host2device(weights);
  else {
    // weights is a 2d array
    // same weight for each subap
    float *tmp;///< Input data
    cutilSafeCall(cudaMalloc((void**)&tmp, sizeof(float)*this->npix*this->npix));
    cutilSafeCall( cudaMemcpy(tmp, weights,sizeof(float)*this->npix*this->npix, cudaMemcpyHostToDevice) );    
    fillweights(this->d_weights->getData(),tmp,this->npix,this->d_weights->getNbElem(),this->device);
    cutilSafeCall( cudaFree(tmp) );
  }

  return EXIT_SUCCESS;
}

int yoga_centroider::init_corr(yoga_wfs *wfs, int isizex, int isizey, float *interpmat)
{
  current_context->set_activeDevice(device);
  if (this->d_corrfnct != 0L) delete this->d_corrfnct;
  if (this->d_corrspot != 0L) delete this->d_corrspot;
  if (this->d_corrnorm != 0L) delete this->d_corrnorm;
  if (this->d_corrmax != 0L) delete this->d_corrmax;
  if (this->d_corr != 0L) delete this->d_corr;
  if (this->d_interpmat != 0L) delete this->d_interpmat;

  this->npix = wfs->npix;

  long *dims_data3 = new long[4];
  dims_data3[0] = 3; 
  dims_data3[1] = 2 * this->npix; dims_data3[2] = 2 * this->npix; dims_data3[3] = this->nvalid; 

  this->d_corrfnct = new yoga_obj<cuFloatComplex>(current_context,dims_data3);
  this->d_corrspot = new yoga_obj<cuFloatComplex>(current_context,dims_data3);

  int mdims[2];
  mdims[0] = (int)dims_data3[1];
  mdims[1] = (int)dims_data3[2];
  cufftHandle *plan=this->d_corrfnct->getPlan();///< FFT plan
  cufftSafeCall(cufftPlanMany(plan, 2 ,mdims,NULL,1,0,NULL,1,0,CUFFT_C2C ,
				(int)dims_data3[3]));

  dims_data3[1] = 2 * this->npix - 1; dims_data3[2] = 2 * this->npix - 1; 
  this->d_corr = new yoga_obj<float>(current_context,dims_data3);

  long *dims_data2 = new long[3];
  dims_data2[0] = 2; 
  dims_data2[1] = 2 * this->npix - 1; dims_data2[2] = 2 * this->npix - 1;
  this->d_corrnorm = new yoga_obj<float>(current_context,dims_data2);

  long *dims_data1 = new long[2];
  dims_data1[0] = 1; 
  dims_data1[1] = this->nvalid;
  this->d_corrmax = new yoga_obj<int>(current_context,dims_data1);

  this->interp_sizex = isizex;
  this->interp_sizey = isizey;
  dims_data2[1] = isizex * isizey; dims_data2[2] = 6;
  this->d_interpmat = new yoga_obj<float>(current_context,dims_data2);
  this->d_interpmat->host2device(interpmat);

  delete[] dims_data1;
  delete[] dims_data2;
  delete[] dims_data3;

  return EXIT_SUCCESS;
}

int yoga_centroider::load_corr(float *corr, float *corr_norm, int ndim)
{
  int nval = (ndim == 3) ? 1 : this->nvalid;

  this->d_corrnorm->host2device(corr_norm);

  float *tmp;///< Input data

  if (ndim == 3) {
    cutilSafeCall(cudaMalloc((void**)&tmp, sizeof(float)*this->npix*this->npix*this->nvalid));
    cutilSafeCall( cudaMemcpy(tmp, corr,sizeof(float)*this->npix*this->npix*this->nvalid, cudaMemcpyHostToDevice) );    
  } else {
    cutilSafeCall(cudaMalloc((void**)&tmp, sizeof(float)*this->npix*this->npix));
    cutilSafeCall( cudaMemcpy(tmp, corr,sizeof(float)*this->npix*this->npix, cudaMemcpyHostToDevice) );    
  }

  fillcorr(this->d_corrfnct->getData(),tmp,this->npix,2*this->npix,this->npix*this->npix*this->nvalid,nval,this->device);

  cutilSafeCall( cudaFree(tmp) );

  yoga_fft(this->d_corrfnct->getData(),this->d_corrfnct->getData(),1,*this->d_corrfnct->getPlan());

  return EXIT_SUCCESS;
}

int yoga_centroider::get_cog(float *cube,float *subsum, float *centroids, int nvalid, int npix, int ntot)
{
  // simple cog

  subap_reduce(ntot,npix*npix,nvalid,cube,subsum);
  
  get_centroids(ntot,npix*npix,nvalid,npix,cube,centroids,subsum,this->scale,this->offset);

  return EXIT_SUCCESS;
}

int yoga_centroider::get_cog(yoga_wfs *wfs, yoga_obj<float> *slopes)
{
  return this->get_cog(wfs->d_bincube->getData(),wfs->d_subsum->getData(),slopes->getData(),
		       wfs->nvalid,wfs->npix,wfs->d_bincube->getNbElem());
}

int yoga_centroider::get_cog(yoga_wfs *wfs)
{
  return this->get_cog(wfs,wfs->d_slopes);
}

int yoga_centroider::get_cog_async(yoga_streams *streams, float *cube,float *subsum, float *centroids, int nvalid, int npix)
{
  // simple cog
  subap_reduce_async(npix*npix,nvalid,streams,cube,subsum);
  get_centroids_async(npix*npix,nvalid,npix,streams,cube,centroids,subsum,this->scale,this->offset);
  return EXIT_SUCCESS;
}

int yoga_centroider::get_cog_async(yoga_wfs *wfs, yoga_obj<float> *slopes)
{
  return this->get_cog_async(wfs->streams,wfs->d_bincube->getData(),wfs->d_subsum->getData(),slopes->getData(),
		       wfs->nvalid,wfs->npix);
}

int yoga_centroider::get_cog_async(yoga_wfs *wfs)
{
  return this->get_cog_async(wfs,wfs->d_slopes);
}

int yoga_centroider::get_tcog(float threshold,float *cube,float *subsum, float *centroids, int nvalid, int npix, int ntot)
{
  // thresholded cog
  // in this case pixels < threshold are put to 0
  subap_reduce(ntot,npix*npix,nvalid,cube,subsum,threshold);
  
  get_centroids(ntot,npix * npix,nvalid,npix,cube,centroids,subsum,threshold,this->scale,this->offset);

  return EXIT_SUCCESS;
}

int yoga_centroider::get_bpcog(int nmax, int npix, int nvalid, float *cube, float *centroids)
{
  // brightest pixels cog
  subap_centromax(npix * npix,nvalid,cube,centroids,npix,nmax,this->scale,this->offset);
  return EXIT_SUCCESS;
}

int yoga_centroider::get_wcog(float *weights,float *cube,float *subsum,float *centroids, int nvalid, int npix, int ntot)
{
  // wcog
  subap_reduce(ntot,npix*npix,nvalid,cube,subsum,weights);
  
  get_centroids(ntot,npix * npix,nvalid,npix,cube,centroids,subsum,weights,this->scale,this->offset);

  return EXIT_SUCCESS;
}


int yoga_centroider::get_tcog(float threshold,yoga_wfs *wfs, float *slopes)
{
  return this->get_tcog(threshold,wfs->d_bincube->getData(),wfs->d_subsum->getData(),
			slopes,wfs->nvalid,wfs->npix,wfs->d_bincube->getNbElem());
}

int yoga_centroider::get_tcog(float threshold,yoga_wfs *wfs)
{
  return this->get_tcog(threshold,wfs,wfs->d_slopes->getData());
}

int yoga_centroider::get_bpcog(int nmax,yoga_wfs *wfs, float *slopes)
{
  return this->get_bpcog(nmax,wfs->npix,wfs->nvalid,wfs->d_bincube->getData(),slopes);
}

int yoga_centroider::get_bpcog(int nmax,yoga_wfs *wfs)
{
  return this->get_bpcog(nmax,wfs,wfs->d_slopes->getData());
}

int yoga_centroider::get_wcog(yoga_wfs *wfs, float *slopes)
{
  return this->get_wcog(this->d_weights->getData(),wfs->d_bincube->getData(),wfs->d_subsum->getData(),
		 slopes,wfs->nvalid,wfs->npix,wfs->d_bincube->getNbElem());
}

int yoga_centroider::get_wcog(yoga_wfs *wfs)
{
  return this->get_wcog(wfs,wfs->d_slopes->getData());
}

int yoga_centroider::get_corr(yoga_wfs *wfs, float *slopes)
{
  //set corrspot to 0
  cutilSafeCall(cudaMemset(this->d_corrspot->getData(), 0,sizeof(cuFloatComplex)*this->d_corrspot->getNbElem()));
  // correlation algorithm
  fillcorr(this->d_corrspot->getData(),wfs->d_bincube->getData(),this->npix,2*this->npix,this->npix*this->npix*this->nvalid,1,this->device);

  yoga_fft(this->d_corrspot->getData(),this->d_corrspot->getData(),1,*this->d_corrfnct->getPlan());

  correl(this->d_corrspot->getData(),this->d_corrfnct->getData(),this->d_corrfnct->getNbElem(),this->device);
  // after this d_corrspot contains the fft of the correl function

  yoga_fft(this->d_corrspot->getData(),this->d_corrspot->getData(),-1,*this->d_corrfnct->getPlan());

  // size is 2 x npix so it is even ...
  roll2real(this->d_corr->getData(),this->d_corrspot->getData(),2*this->npix,(2*this->npix)*(2*this->npix),
	    this->d_corrspot->getNbElem(),this->device);
  //here need to normalize
  corr_norm(this->d_corr->getData(),this->d_corrnorm->getData(),this->d_corrnorm->getNbElem(),this->d_corr->getNbElem(),this->device);

  // need to find max for each subap
  // if the corr array for one subap is greater than 20x20
  // it won't fit in shared mem
  // so we window around the center of the array, 
  // the max is expected to be found inside the npix x npix central part anyway  
  int nbmax = (2*this->npix-1 > 20) ? this->npix : 2*this->npix-1;
  int xoff  = this->d_corr->getDims(1)/2 - nbmax/2 ;
  int yoff  = this->d_corr->getDims(2)/2 - nbmax/2;

  subap_sortmaxi(nbmax * nbmax,this->nvalid,this->d_corr->getData(),this->d_corrmax->getData() ,1,
		 xoff,yoff,nbmax,this->d_corr->getDims(1));

  // do parabolic interpolation
  subap_pinterp(this->interp_sizex * this->interp_sizey,this->nvalid,this->d_corr->getData(),
		this->d_corrmax->getData(),slopes,this->d_interpmat->getData(),
		this->interp_sizex ,this->interp_sizey,this->nvalid,2*this->npix-1,this->scale,this->offset);

  return EXIT_SUCCESS;
}

int yoga_centroider::get_corr(yoga_wfs *wfs)
{
  return this->get_corr(wfs,wfs->d_slopes->getData());

}

// deprecated
int yoga_centroider::init_nmax(int nmax)
{
  if (this->d_validpix != 0L) delete this->d_validpix;
  if (this->d_validindx!= 0L)  delete  this->d_validindx;

  long *dims_data2 = new long[3];
  dims_data2[0] = 2; 
  dims_data2[1] = nmax; dims_data2[2] = this->nvalid; 

  this->d_validpix = new yoga_obj<float>(current_context,dims_data2);
  this->d_validindx = new yoga_obj<int>(current_context,dims_data2);

  delete[] dims_data2;

  return EXIT_SUCCESS;
}

int yoga_centroider::get_pyr(float *cube,float *subsum, float *centroids, int *subindx, int *subindy, int nvalid, int ns, int nim)
{
  // simple cog
  pyr_slopes(centroids,cube,subindx,subindy,subsum,ns,nvalid,nim,this->device);
  
  return EXIT_SUCCESS;
}

int yoga_centroider::get_pyr(yoga_wfs *wfs, yoga_obj<float> *slopes)
{
  return this->get_pyr(wfs->d_bincube->getData(),wfs->d_subsum->getData(),slopes->getData(),
		       wfs->d_validsubsx->getData(),wfs->d_validsubsy->getData(), wfs->nvalid,
		       wfs->nfft/wfs->nrebin,4);
}

int yoga_centroider::get_pyr(yoga_wfs *wfs)
{
  return this->get_pyr(wfs,wfs->d_slopes);
}
