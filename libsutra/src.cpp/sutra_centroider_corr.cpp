#include <sutra_centroider_corr.h>
#include <string>

sutra_centroider_corr::sutra_centroider_corr(carma_context *context,
                                             sutra_sensors *sensors, int nwfs,
                                             long nvalid, float offset,
                                             float scale, int device)
    : sutra_centroider(context, sensors, nwfs, nvalid, offset, scale, device) {
  context->set_activeDevice(device, 1);

  this->d_corrfnct = 0L;
  this->d_corrspot = 0L;
  this->d_corrnorm = 0L;
  this->d_corrmax = 0L;
  this->d_corr = 0L;
  this->d_interpmat = 0L;

  this->npix = 0;
  this->interp_sizex = 0;
  this->interp_sizey = 0;
}

sutra_centroider_corr::~sutra_centroider_corr() {
  if (this->d_corrfnct != 0L) delete this->d_corrfnct;
  if (this->d_corrspot != 0L) delete this->d_corrspot;
  if (this->d_corrnorm != 0L) delete this->d_corrnorm;
  if (this->d_corrmax != 0L) delete this->d_corrmax;
  if (this->d_corr != 0L) delete this->d_corr;
  if (this->d_interpmat != 0L) delete this->d_interpmat;
}

string sutra_centroider_corr::get_type() { return "corr"; }

int sutra_centroider_corr::init_bincube(int npix) {
  this->npix = npix;

  return EXIT_SUCCESS;
}

int sutra_centroider_corr::init_corr(int isizex, int isizey, float *interpmat) {
  current_context->set_activeDevice(device, 1);
  if (this->d_corrfnct != 0L) delete this->d_corrfnct;
  if (this->d_corrspot != 0L) delete this->d_corrspot;
  if (this->d_corrnorm != 0L) delete this->d_corrnorm;
  if (this->d_corrmax != 0L) delete this->d_corrmax;
  if (this->d_corr != 0L) delete this->d_corr;
  if (this->d_interpmat != 0L) delete this->d_interpmat;

  this->npix = wfs->npix;

  long *dims_data3 = new long[4];
  dims_data3[0] = 3;
  dims_data3[1] = 2 * this->npix;
  dims_data3[2] = 2 * this->npix;
  dims_data3[3] = this->nvalid;

  this->d_corrfnct = new carma_obj<cuFloatComplex>(current_context, dims_data3);
  this->d_corrspot = new carma_obj<cuFloatComplex>(current_context, dims_data3);

  int mdims[2];
  mdims[0] = (int)dims_data3[1];
  mdims[1] = (int)dims_data3[2];
  cufftHandle *plan = this->d_corrfnct->getPlan();  ///< FFT plan
  carmafftSafeCall(cufftPlanMany(plan, 2, mdims, NULL, 1, 0, NULL, 1, 0,
                                 CUFFT_C2C, (int)dims_data3[3]));

  dims_data3[1] = 2 * this->npix - 1;
  dims_data3[2] = 2 * this->npix - 1;
  this->d_corr = new carma_obj<float>(current_context, dims_data3);

  long *dims_data2 = new long[3];
  dims_data2[0] = 2;
  dims_data2[1] = 2 * this->npix - 1;
  dims_data2[2] = 2 * this->npix - 1;
  this->d_corrnorm = new carma_obj<float>(current_context, dims_data2);

  long *dims_data1 = new long[2];
  dims_data1[0] = 1;
  dims_data1[1] = this->nvalid;
  this->d_corrmax = new carma_obj<int>(current_context, dims_data1);

  this->interp_sizex = isizex;
  this->interp_sizey = isizey;
  dims_data2[1] = isizex * isizey;
  dims_data2[2] = 6;
  this->d_interpmat = new carma_obj<float>(current_context, dims_data2);
  this->d_interpmat->host2device(interpmat);

  delete[] dims_data1;
  delete[] dims_data2;
  delete[] dims_data3;

  return EXIT_SUCCESS;
}

int sutra_centroider_corr::load_corr(float *corr, float *corr_norm, int ndim) {
  current_context->set_activeDevice(device, 1);
  int nval = (ndim == 3) ? 1 : this->nvalid;

  this->d_corrnorm->host2device(corr_norm);

  float *tmp;  ///< Input data

  if (ndim == 3) {
    carmaSafeCall(cudaMalloc(
        (void **)&tmp, sizeof(float) * this->npix * this->npix * this->nvalid));
    carmaSafeCall(cudaMemcpy(
        tmp, corr, sizeof(float) * this->npix * this->npix * this->nvalid,
        cudaMemcpyHostToDevice));
  } else {
    carmaSafeCall(
        cudaMalloc((void **)&tmp, sizeof(float) * this->npix * this->npix));
    carmaSafeCall(cudaMemcpy(tmp, corr, sizeof(float) * this->npix * this->npix,
                             cudaMemcpyHostToDevice));
  }

  fillcorr(*(this->d_corrfnct), tmp, this->npix, 2 * this->npix,
           this->npix * this->npix * this->nvalid, nval,
           this->current_context->get_device(device));

  carmaSafeCall(cudaFree(tmp));

  carma_fft<cuFloatComplex, cuFloatComplex>(*(this->d_corrfnct),
                                            *(this->d_corrfnct), 1,
                                            *this->d_corrfnct->getPlan());

  return EXIT_SUCCESS;
}

int sutra_centroider_corr::get_cog(carma_streams *streams, float *cube,
                                   float *subsum, float *centroids, int nvalid,
                                   int npix, int ntot) {
  current_context->set_activeDevice(device, 1);
  // set corrspot to 0
  carmaSafeCall(
      cudaMemset(*(this->d_corrspot), 0,
                 sizeof(cuFloatComplex) * this->d_corrspot->getNbElem()));
  // correlation algorithm

  fillcorr(*(this->d_corrspot), cube, this->npix, 2 * this->npix,
           this->npix * this->npix * this->nvalid, 1,
           this->current_context->get_device(device));

  carma_fft<cuFloatComplex, cuFloatComplex>(*(this->d_corrspot),
                                            *(this->d_corrspot), 1,
                                            *this->d_corrfnct->getPlan());

  correl(*(this->d_corrspot), *(this->d_corrfnct),
         this->d_corrfnct->getNbElem(),
         this->current_context->get_device(device));
  // after this d_corrspot contains the fft of the correl function

  carma_fft<cuFloatComplex, cuFloatComplex>(*(this->d_corrspot),
                                            *(this->d_corrspot), -1,
                                            *this->d_corrfnct->getPlan());

  // size is 2 x npix so it is even ...
  roll2real(*(this->d_corr), *(this->d_corrspot), 2 * this->npix,
            (2 * this->npix) * (2 * this->npix), this->d_corrspot->getNbElem(),
            this->current_context->get_device(device));
  // here need to normalize
  corr_norm(*(this->d_corr), *(this->d_corrnorm), this->d_corrnorm->getNbElem(),
            this->d_corr->getNbElem(),
            this->current_context->get_device(device));

  // need to find max for each subap
  // if the corr array for one subap is greater than 20x20
  // it won't fit in shared mem
  // so we window around the center of the array,
  // the max is expected to be found inside the npix x npix central part anyway
  int nbmax = (2 * this->npix - 1 > 20) ? this->npix : 2 * this->npix - 1;
  int xoff = this->d_corr->getDims(1) / 2 - nbmax / 2;
  int yoff = this->d_corr->getDims(2) / 2 - nbmax / 2;

  subap_sortmaxi<float>(nbmax * nbmax, this->nvalid, *(this->d_corr),
                        *(this->d_corrmax), 1, xoff, yoff, nbmax,
                        this->d_corr->getDims(1));

  // do parabolic interpolation
  subap_pinterp<float>(this->interp_sizex * this->interp_sizey, this->nvalid,
                       *(this->d_corr), *(this->d_corrmax), centroids,
                       *(this->d_interpmat), this->interp_sizex,
                       this->interp_sizey, this->nvalid, 2 * this->npix - 1,
                       this->scale, this->offset);
  return EXIT_SUCCESS;
}

int sutra_centroider_corr::get_cog(float *subsum, float *slopes, bool noise) {
  if (this->wfs != nullptr) {
    return this->get_cog(wfs->streams, *wfs->d_bincube, subsum, slopes,
                         wfs->nvalid, wfs->npix, wfs->d_bincube->getNbElem());
  }

  DEBUG_TRACE("this->wfs was not initialized");
  return EXIT_FAILURE;
}

int sutra_centroider_corr::get_cog() {
  if (this->wfs != nullptr)
    return this->get_cog(*(wfs->d_subsum), *(wfs->d_slopes), true);

  DEBUG_TRACE("this->wfs was not initialized");
  return EXIT_FAILURE;
}
