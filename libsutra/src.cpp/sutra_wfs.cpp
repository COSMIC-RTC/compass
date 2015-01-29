#include <sutra_wfs.h>
#include <sutra_ao_utils.h>
#include <carma_utils.h>

sutra_wfs::sutra_wfs(carma_context *context, const char* type, long nxsub,
    long nvalid, long npix, long nphase, long nrebin, long nfft, long ntot,
    long npup, float pdiam, float nphotons, int lgs, int device) {
  this->d_camplipup = 0L;
  this->d_camplifoc = 0L;
  this->d_fttotim = 0L;
  this->d_ftkernel = 0L;
  this->d_pupil = 0L;
  this->d_hrimg = 0L;
  this->d_bincube = 0L;
  this->d_binimg = 0L;
  this->d_subsum = 0L;
  this->d_offsets = 0L;
  this->d_fluxPerSub = 0L;
  this->d_sincar = 0L;
  this->d_submask = 0L;
  this->d_hrmap = 0L;
  this->d_isvalid = 0L;
  this->d_slopes = 0L;
  this->image_telemetry = 0L;
  this->d_phasemap = 0L;
  this->d_binmap = 0L;
  this->d_validsubsx = 0L;
  this->d_validsubsy = 0L;
  this->d_istart = 0L;
  this->d_jstart = 0L;
  this->d_psum = 0L;
  this->d_phalfxy = 0L;
  this->d_poffsets = 0L;
  this->pyr_cx = 0L;
  this->pyr_cy = 0L;

  this->d_gs = 0L;

  this->current_context = context;

  this->type = type;

  this->noise = 0;
  this->nxsub = nxsub;
  this->nvalid = nvalid;
  this->npix = npix;
  this->nphase = nphase;
  this->nrebin = nrebin;
  this->nfft = nfft;
  this->ntot = ntot;
  this->npup = npup;
  this->subapd = pdiam;
  this->nphot = nphotons;
  this->lgs = (lgs == 1 ? true : false);
  this->device = device;
  context->set_activeDevice(device);
  this->nmaxhr = nvalid;
  this->nffthr = 1;

  this->kernconv = false;

  long *dims_data1 = new long[2];
  dims_data1[0] = 1;
  long *dims_data2 = new long[3];
  dims_data2[0] = 2;
  long *dims_data3 = new long[4];
  dims_data3[0] = 3;

  dims_data2[1] = npix * nxsub;
  dims_data2[2] = npix * nxsub;
  if ((this->type == "pyr") || (this->type == "pyr3")) {
    dims_data2[1] += 3;
    dims_data2[2] += 3;
  }

  this->d_binimg = new carma_obj<float>(context, dims_data2);
  // using 1 stream for telemetry
  this->image_telemetry = new carma_host_obj<float>(dims_data2, MA_PAGELOCK, 1);

  dims_data3[1] = nfft;
  dims_data3[2] = nfft;
  if (this->type == "pyr")
    dims_data3[3] = 4;
  if (this->type == "pyr3")
    dims_data3[3] = 3;

  int nslaves = 0;

  if (this->type == "sh") {
    if (nslaves > 0)
      dims_data3[3] = nvalid / (nslaves + 1);
    else
      dims_data3[3] = nvalid;
  }
  else
	  this->d_hrimg = new carma_obj<float>(context, dims_data3); // Useless for SH

  int mdims[2];
  if (this->type == "sh") {
    //this->d_submask = new carma_obj<float>(context, dims_data3); // Useless for SH
    this->d_camplipup = new carma_obj<cuFloatComplex>(context, dims_data3);
    this->d_camplifoc = new carma_obj<cuFloatComplex>(context, dims_data3);
    mdims[0] = (int) dims_data3[1];
    mdims[1] = (int) dims_data3[2];
    cufftHandle *plan = this->d_camplipup->getPlan(); ///< FFT plan
    cufftSafeCall(
        cufftPlanMany(plan, 2 ,mdims,NULL,1,0,NULL,1,0, CUFFT_C2C ,(int)dims_data3[3]));
    dims_data3[1] = npix;
    dims_data3[2] = npix;
  }

  if ((this->type == "pyr") || (this->type == "pyr3")) {
    dims_data2[1] = ntot;
    dims_data2[2] = ntot;
    this->d_submask = new carma_obj<float>(context, dims_data2);
    this->d_camplipup = new carma_obj<cuFloatComplex>(context, dims_data2);
    this->d_camplifoc = new carma_obj<cuFloatComplex>(context, dims_data2);
    cufftHandle *plan = this->d_camplipup->getPlan(); ///< FFT plan
    cufftSafeCall(cufftPlan2d(plan, dims_data2[1], dims_data2[2], CUFFT_C2C));

    dims_data3[1] = nfft / nrebin;
    dims_data3[2] = nfft / nrebin;
  }

  this->d_bincube = new carma_obj<float>(context, dims_data3);

  this->nstreams = 1;
  while (nvalid % this->nstreams != 0)
    nstreams--;
  cerr << "wfs uses " << nstreams << " streams" << endl;
  this->streams = new carma_streams(nstreams);

  dims_data1[1] = 2 * nvalid;
  this->d_slopes = new carma_obj<float>(context, dims_data1);

  if (this->type == "sh") {
    if (this->ntot != this->nfft) {
      // this is the big array => we use nmaxhr and treat it sequentially
      int mnmax = 500;
      if (nvalid > 2 * mnmax) {
        this->nmaxhr = mnmax;
        int tmp0 = nvalid % mnmax;
        int tmp = 0;
        for (int cc = 1; cc < mnmax / 5; cc++) {
          tmp = nvalid % (mnmax + cc);
          if ((tmp > tmp0) || (tmp == 0)) {
            if (tmp == 0)
              tmp0 = 2 * mnmax;
            else
              tmp = tmp0;
            this->nmaxhr = mnmax + cc;
          }
        }
        this->nffthr = (
            nvalid % this->nmaxhr == 0 ?
                nvalid / this->nmaxhr : nvalid / this->nmaxhr + 1);
      }

      dims_data3[1] = ntot;
      dims_data3[2] = ntot;
      dims_data3[3] = nmaxhr;
      this->d_fttotim = new carma_obj<cuFloatComplex>(context, dims_data3);
      mdims[0] = (int) dims_data3[1];
      mdims[1] = (int) dims_data3[2];
      cufftHandle *plan = this->d_fttotim->getPlan(); ///< FFT plan
      cufftSafeCall(
          cufftPlanMany(plan, 2 ,mdims,NULL,1,0,NULL,1,0,CUFFT_C2C , (int)dims_data3[3]));

      dims_data1[1] = nfft * nfft;
      this->d_hrmap = new carma_obj<int>(context, dims_data1);

    } else {
      if (this->lgs) {
        dims_data3[1] = ntot;
        dims_data3[2] = ntot;
        dims_data3[3] = nvalid;
        this->d_fttotim = new carma_obj<cuFloatComplex>(context, dims_data3);
        mdims[0] = (int) dims_data3[1];
        mdims[1] = (int) dims_data3[2];
        cufftHandle *plan = this->d_fttotim->getPlan(); ///< FFT plan
        cufftSafeCall(
            cufftPlanMany(plan, 2 ,mdims,NULL,1,0,NULL,1,0,CUFFT_C2C , (int)dims_data3[3]));
      }
    }

    dims_data2[1] = ntot;
    dims_data2[2] = ntot;
    this->d_ftkernel = new carma_obj<cuFloatComplex>(context, dims_data2);

    dims_data2[1] = npup;
    dims_data2[2] = npup;
    this->d_pupil = new carma_obj<float>(context, dims_data2);

    dims_data2[1] = nphase;
    dims_data2[2] = nphase;
    this->d_offsets = new carma_obj<float>(context, dims_data2);

    dims_data1[1] = nxsub;
    this->d_istart = new carma_obj<int>(context, dims_data1);
    this->d_jstart = new carma_obj<int>(context, dims_data1);

    dims_data2[1] = nrebin * nrebin;
    dims_data2[2] = npix * npix;
    this->d_binmap = new carma_obj<int>(context, dims_data2);
  }

  if (this->type == "pyr") {
    dims_data3[1] = nfft;
    dims_data3[2] = nfft;
    dims_data3[3] = 4;
    this->d_fttotim = new carma_obj<cuFloatComplex>(context, dims_data3);
    mdims[0] = (int) dims_data3[1];
    mdims[1] = (int) dims_data3[2];
    cufftHandle *plan = this->d_fttotim->getPlan(); ///< FFT plan
    cufftSafeCall(
        cufftPlanMany(plan, 2 ,mdims,NULL,1,0,NULL,1,0,CUFFT_C2C , (int)dims_data3[3]));

    dims_data2[1] = ntot;
    dims_data2[2] = ntot;
    this->d_pupil = new carma_obj<float>(context, dims_data2);

    dims_data1[1] = npup;
    this->pyr_cx = new carma_host_obj<int>(dims_data1, MA_WRICOMB);
    this->pyr_cy = new carma_host_obj<int>(dims_data1, MA_WRICOMB);
    dims_data2[1] = nfft;
    dims_data2[2] = nfft;
    this->d_poffsets = new carma_obj<cuFloatComplex>(context, dims_data2);
    dims_data2[1] = ntot;
    dims_data2[2] = ntot;
    this->d_phalfxy = new carma_obj<cuFloatComplex>(context, dims_data2);
    dims_data2[1] = nfft;
    dims_data2[2] = nfft;
    this->d_sincar = new carma_obj<float>(context, dims_data2);
  }

  dims_data1[1] = nvalid;
  this->d_subsum = new carma_obj<float>(context, dims_data1);

  if (this->type == "pyr")
    this->d_psum = new carma_obj<float>(context, dims_data1);

  this->d_fluxPerSub = new carma_obj<float>(context, dims_data1);
  this->d_validsubsx = new carma_obj<int>(context, dims_data1);
  this->d_validsubsy = new carma_obj<int>(context, dims_data1);

  dims_data2[1] = nxsub;
  dims_data2[2] = nxsub;
  this->d_isvalid = new carma_obj<int>(context, dims_data2);

  dims_data2[1] = nphase * nphase;
  dims_data2[2] = nvalid;
  this->d_phasemap = new carma_obj<int>(context, dims_data2);
}

sutra_wfs::sutra_wfs(carma_context *context, long nxsub, long nvalid,
    long nphase, long npup, float pdiam, int device) {
  this->d_camplipup = 0L;
  this->d_camplifoc = 0L;
  this->d_fttotim = 0L;
  this->d_ftkernel = 0L;
  this->d_pupil = 0L;
  this->d_hrimg = 0L;
  this->d_bincube = 0L;
  this->d_binimg = 0L;
  this->d_subsum = 0L;
  this->d_offsets = 0L;
  this->d_fluxPerSub = 0L;
  this->d_sincar = 0L;
  this->d_submask = 0L;
  this->d_hrmap = 0L;
  this->d_isvalid = 0L;
  this->d_slopes = 0L;
  this->image_telemetry = 0L;
  this->d_phasemap = 0L;
  this->d_binmap = 0L;
  this->d_validsubsx = 0L;
  this->d_validsubsy = 0L;
  this->d_istart = 0L;
  this->d_jstart = 0L;
  this->d_psum = 0L;
  this->d_phalfxy = 0L;
  this->d_poffsets = 0L;
  this->pyr_cx = 0L;
  this->pyr_cy = 0L;

  this->current_context = context;

  this->type = "geo";

  this->kernconv = false;
  this->d_gs = 0L;
  this->noise = 0;
  this->nxsub = nxsub;
  this->nvalid = nvalid;
  this->nphase = nphase;
  this->npup = npup;
  this->subapd = pdiam;
  this->device = device;
  context->set_activeDevice(device);

  this->npix = 0;
  this->nrebin = 0;
  this->nfft = 0;
  this->ntot = 0;
  this->nphot = 0;
  this->lgs = false;
  this->nmaxhr = 0;
  this->nffthr = 0;

  this->nstreams = 1; //nvalid/10;
  this->streams = new carma_streams(nstreams);

  long *dims_data1 = new long[2];
  dims_data1[0] = 1;
  long *dims_data2 = new long[3];
  dims_data2[0] = 2;
  long *dims_data3 = new long[4];
  dims_data3[0] = 3;

  dims_data1[1] = 2 * nvalid;
  this->d_slopes = new carma_obj<float>(context, dims_data1);

  dims_data2[1] = npup;
  dims_data2[2] = npup;
  this->d_pupil = new carma_obj<float>(context, dims_data2);

  dims_data2[1] = nphase;
  dims_data2[2] = nphase;
  this->d_offsets = new carma_obj<float>(context, dims_data2);

  dims_data1[1] = nxsub;
  this->d_istart = new carma_obj<int>(context, dims_data1);
  this->d_jstart = new carma_obj<int>(context, dims_data1);

  dims_data1[1] = nvalid;
  this->d_subsum = new carma_obj<float>(context, dims_data1);

  this->d_fluxPerSub = new carma_obj<float>(context, dims_data1);
  this->d_validsubsx = new carma_obj<int>(context, dims_data1);
  this->d_validsubsy = new carma_obj<int>(context, dims_data1);

  dims_data2[1] = nxsub;
  dims_data2[2] = nxsub;
  this->d_isvalid = new carma_obj<int>(context, dims_data2);

  dims_data2[1] = nphase * nphase;
  dims_data2[2] = nvalid;
  this->d_phasemap = new carma_obj<int>(context, dims_data2);
}

sutra_wfs::~sutra_wfs() {
  current_context->set_activeDevice(device);
  if (this->d_camplipup != 0L)
    delete this->d_camplipup;
  if (this->d_camplifoc != 0L)
    delete this->d_camplifoc;

  if (this->d_fttotim != 0L)
    delete this->d_fttotim;

  if (this->d_ftkernel != 0L)
    delete this->d_ftkernel;

  if (this->d_pupil != 0L)
    delete this->d_pupil;
  if (this->d_hrimg != 0L)
    delete this->d_hrimg;
  if (this->d_bincube != 0L)
    delete this->d_bincube;
  if (this->d_binimg != 0L)
    delete this->d_binimg;
  if (this->d_subsum != 0L)
    delete this->d_subsum;
  if (this->d_offsets != 0L)
    delete this->d_offsets;
  if (this->d_fluxPerSub != 0L)
    delete this->d_fluxPerSub;
  if (this->d_sincar != 0L)
    delete this->d_sincar;
  if (this->d_submask != 0L)
    delete this->d_submask;
  if (this->d_hrmap != 0L)
    delete this->d_hrmap;

  if (this->d_isvalid != 0L)
    delete this->d_isvalid;
  if (this->d_slopes != 0L)
    delete this->d_slopes;

  if (this->image_telemetry != 0L)
    delete this->image_telemetry;

  if (this->d_phasemap != 0L)
    delete this->d_phasemap;
  if (this->d_binmap != 0L)
    delete this->d_binmap;
  if (this->d_validsubsx != 0L)
    delete this->d_validsubsx;
  if (this->d_validsubsy != 0L)
    delete this->d_validsubsy;
  if (this->d_istart != 0L)
    delete this->d_istart;
  if (this->d_jstart != 0L)
    delete this->d_jstart;

  if (this->d_psum != 0L)
    delete this->d_psum;
  if (this->d_phalfxy != 0L)
    delete this->d_phalfxy;
  if (this->d_poffsets != 0L)
    delete this->d_poffsets;
  if (this->pyr_cx != 0L)
    delete this->pyr_cx;
  if (this->pyr_cy != 0L)
    delete this->pyr_cy;

  if (this->lgs)
    delete this->d_gs->d_lgs;

  delete this->d_gs;

  delete this->streams;

  //delete this->current_context;
}

int sutra_wfs::wfs_initgs(float xpos, float ypos, float lambda, float mag,
    long size, float noise, long seed) {
  this->d_gs = new sutra_source(current_context, xpos, ypos, lambda, mag, size,
      "wfs", this->device);
  this->noise = noise;
  if (noise > -1) {
    this->d_bincube->init_prng(seed);
    this->d_bincube->prng('N', noise, 0.0f);
  }
  if (noise > 0) {
    this->d_binimg->init_prng(seed);
    this->d_binimg->prng('N', noise, 0.0f);
  }

  if (this->lgs) {
    this->d_gs->d_lgs = new sutra_lgs(current_context, this->nvalid, this->ntot,
        this->nmaxhr);
    this->d_gs->lgs = this->lgs;
  }

  return EXIT_SUCCESS;
}

int sutra_wfs::wfs_initarrays(int *phasemap, int *hrmap, int *binmap,
    float *offsets, float *pupil, float *fluxPerSub, int *isvalid,
    int *validsubsx, int *validsubsy, int *istart, int *jstart,
    cuFloatComplex *kernel) {
  this->d_phasemap->host2device(phasemap);
  this->d_offsets->host2device(offsets);
  this->d_pupil->host2device(pupil);
  this->d_binmap->host2device(binmap);
  this->d_fluxPerSub->host2device(fluxPerSub);
  if (this->ntot != this->nfft)
    this->d_hrmap->host2device(hrmap);
  this->d_validsubsx->host2device(validsubsx);
  this->d_validsubsy->host2device(validsubsy);
  this->d_isvalid->host2device(isvalid);
  this->d_istart->host2device(istart);
  this->d_jstart->host2device(jstart);
  this->d_ftkernel->host2device(kernel);

  return EXIT_SUCCESS;
}

int sutra_wfs::wfs_initarrays(int *phasemap, float *offsets, float *pupil,
    float *fluxPerSub, int *isvalid, int *validsubsx, int *validsubsy,
    int *istart, int *jstart) {
  this->d_phasemap->host2device(phasemap);
  this->d_offsets->host2device(offsets);
  this->d_pupil->host2device(pupil);
  this->d_fluxPerSub->host2device(fluxPerSub);
  this->d_validsubsx->host2device(validsubsx);
  this->d_validsubsy->host2device(validsubsy);
  this->d_isvalid->host2device(isvalid);
  this->d_istart->host2device(istart);
  this->d_jstart->host2device(jstart);

  return EXIT_SUCCESS;
}

int sutra_wfs::wfs_initarrays(cuFloatComplex *halfxy, cuFloatComplex *offsets,
    float *focmask, float *pupil, int *isvalid, int *cx, int *cy, float *sincar,
    int *phasemap, int *validsubsx, int *validsubsy) {
  this->d_phalfxy->host2device(halfxy);
  this->d_poffsets->host2device(offsets);
  if (this->type != "sh")
	  this->d_submask->host2device(focmask);
  this->d_pupil->host2device(pupil);
  this->pyr_cx->fill_from(cx);
  this->pyr_cy->fill_from(cy);
  this->d_isvalid->host2device(isvalid);
  this->d_sincar->host2device(sincar);
  this->d_phasemap->host2device(phasemap);
  this->d_validsubsx->host2device(validsubsx);
  this->d_validsubsy->host2device(validsubsy);

  return EXIT_SUCCESS;
}

int sutra_wfs::load_kernels(float *lgskern) {
  if (this->lgs)
    this->d_gs->d_lgs->load_kernels(lgskern, this->current_context->get_device(device));

  return EXIT_SUCCESS;
}

int sutra_wfs::sensor_trace(sutra_atmos *yatmos) {
  //do raytracing to get the phase
  this->d_gs->raytrace(yatmos);

  // same with textures
  // apparently not working for large sizes...
  //this->d_gs->raytrace_shm(yatmos);

  return EXIT_SUCCESS;
}

int sutra_wfs::sensor_trace(sutra_dms *ydm, int rst) {
  //do raytracing to get the phase
  this->d_gs->raytrace(ydm, rst);

  return EXIT_SUCCESS;
}

/////////////////////////////////////////////////////////
// COMPUTATION OF THE SHACK-HARTMANN WAVEFRONT SENSOR  //
/////////////////////////////////////////////////////////

int sutra_wfs::comp_sh_generic() {
  current_context->set_activeDevice(device);

  // segment phase and fill cube of complex ampli with exp(i*phase_seg)
  fillcamplipup(this->d_camplipup->getData(),
      this->d_gs->d_phase->d_screen->getData(), this->d_offsets->getData(),
      this->d_pupil->getData(), this->d_gs->scale, this->d_istart->getData(),
      this->d_jstart->getData(), this->d_validsubsx->getData(),
      this->d_validsubsy->getData(), this->nphase,
      this->d_gs->d_phase->d_screen->getDims(1), this->nfft,
      this->nphase * this->nphase * this->nvalid, this->current_context->get_device(device));

  // do fft of the cube  
  carma_fft(this->d_camplipup->getData(), this->d_camplifoc->getData(), 1,
      *this->d_camplipup->getPlan());
  // get the hrimage by taking the | |^2
  // keep it in amplifoc to save mem space
  abs2c(this->d_camplifoc->getData(), this->d_camplifoc->getData(),
      this->d_camplifoc->getDims(1) * this->d_camplifoc->getDims(2)
          * this->d_camplifoc->getDims(3), current_context->get_device(device));

  //set bincube to 0 or noise
  cutilSafeCall(
      cudaMemset(this->d_bincube->getData(), 0,
          sizeof(float) * this->d_bincube->getNbElem()));
  // increase fov if required
  // and fill bincube with data from hrimg
  // we need to do this sequentially if nvalid > nmaxhr to
  // keep raesonable mem occupancy
  if (this->ntot != this->nfft) {

    for (int cc = 0; cc < this->nffthr; cc++) {
      cutilSafeCall(
          cudaMemset(this->d_fttotim->getData(), 0,
              sizeof(cuFloatComplex) * this->d_fttotim->getNbElem()));

      int indxstart1, indxstart2=0, indxstart3;

      if ((cc == this->nffthr - 1) && (this->nvalid % this->nmaxhr != 0)) {
        indxstart1 = this->d_camplifoc->getNbElem()
            - this->nfft * this->nfft * this->nmaxhr;
        if (this->lgs)
          indxstart2 = this->ntot * this->nvalid - this->ntot * this->nmaxhr;
        indxstart3 = this->d_bincube->getNbElem()
            - this->npix * this->npix * this->nmaxhr;
      } else {
        indxstart1 = this->nfft * this->nfft * this->nmaxhr * cc;
        if (this->lgs)
          indxstart2 = this->ntot * this->nmaxhr * cc;
        indxstart3 = this->npix * this->npix * this->nmaxhr * cc;
      }

      cuFloatComplex *data = this->d_camplifoc->getData();
      indexfill(this->d_fttotim->getData(), &(data[indxstart1]),
          this->d_hrmap->getData(), this->nfft, this->ntot,
          this->nfft * this->nfft * this->nmaxhr, this->current_context->get_device(device));

      if (this->lgs) {
        // compute lgs spot on the fly from binned profile image
        this->d_gs->d_lgs->lgs_makespot(this->current_context->get_device(device), indxstart2);
        // convolve with psf
        carma_fft(this->d_fttotim->getData(), this->d_fttotim->getData(), 1,
            *this->d_fttotim->getPlan());

        convolve(this->d_fttotim->getData(),
            this->d_gs->d_lgs->d_ftlgskern->getData(),
            this->d_fttotim->getNbElem(), this->current_context->get_device(device));

        carma_fft(this->d_fttotim->getData(), this->d_fttotim->getData(), -1,
            *this->d_fttotim->getPlan());

      }

      if (this->kernconv) {
        carma_fft(this->d_fttotim->getData(), this->d_fttotim->getData(), 1,
            *this->d_fttotim->getPlan());

        convolve_cube(this->d_fttotim->getData(), this->d_ftkernel->getData(),
            this->d_fttotim->getNbElem(), this->d_ftkernel->getNbElem(),
            this->current_context->get_device(device));

        carma_fft(this->d_fttotim->getData(), this->d_fttotim->getData(), -1,
            *this->d_fttotim->getPlan());

      }

      float *data2 = this->d_bincube->getData();
      //fprintf(stderr, "[%s@%d]: I'm here!\n", __FILE__, __LINE__);
      if (this->nstreams > 1) {
        fillbincube_async(this->streams, &(data2[indxstart3]),
            this->d_fttotim->getData(), this->d_binmap->getData(),
            this->ntot * this->ntot, this->npix * this->npix,
            this->nrebin * this->nrebin, this->nmaxhr, this->current_context->get_device(device));
      } else {
        fillbincube(&(data2[indxstart3]), this->d_fttotim->getData(),
            this->d_binmap->getData(), this->ntot * this->ntot,
            this->npix * this->npix, this->nrebin * this->nrebin, this->nmaxhr,
            this->current_context->get_device(device));
      //fprintf(stderr, "[%s@%d]: I'm here!\n", __FILE__, __LINE__);
      }
    }
  } else {
    if (this->lgs) {
      this->d_gs->d_lgs->lgs_makespot(this->current_context->get_device(device), 0);

      carma_fft(this->d_camplifoc->getData(), this->d_fttotim->getData(), 1,
          *this->d_fttotim->getPlan());

      convolve(this->d_fttotim->getData(),
          this->d_gs->d_lgs->d_ftlgskern->getData(),
          this->d_fttotim->getNbElem(), this->current_context->get_device(device));

      carma_fft(this->d_fttotim->getData(), this->d_fttotim->getData(), -1,
          *this->d_fttotim->getPlan());

      if (this->nstreams > 1){
        fillbincube_async(this->streams, this->d_bincube->getData(),
            this->d_fttotim->getData(), this->d_binmap->getData(),
            this->nfft * this->nfft, this->npix * this->npix,
            this->nrebin * this->nrebin, this->nvalid, this->current_context->get_device(device));
      } else {
        fillbincube(this->d_bincube->getData(), this->d_fttotim->getData(),
            this->d_binmap->getData(), this->nfft * this->nfft,
            this->npix * this->npix, this->nrebin * this->nrebin, this->nvalid,
            this->current_context->get_device(device));
      }
    } else {
      if (this->kernconv) {
        carma_fft(this->d_camplifoc->getData(), this->d_camplifoc->getData(), 1,
            *this->d_camplipup->getPlan());

        convolve_cube(this->d_camplifoc->getData(), this->d_ftkernel->getData(),
            this->d_camplifoc->getNbElem(), this->d_ftkernel->getNbElem(),
            this->current_context->get_device(device));

        carma_fft(this->d_camplifoc->getData(), this->d_camplifoc->getData(),
            -1, *this->d_camplipup->getPlan());
      }

      if (this->nstreams > 1) {
        fillbincube_async(this->streams, this->d_bincube->getData(),
            this->d_camplifoc->getData(), this->d_binmap->getData(),
            this->nfft * this->nfft, this->npix * this->npix,
            this->nrebin * this->nrebin, this->nvalid, this->current_context->get_device(device));
      } else {
        fillbincube(this->d_bincube->getData(), this->d_camplifoc->getData(),
            this->d_binmap->getData(), this->nfft * this->nfft,
            this->npix * this->npix, this->nrebin * this->nrebin, this->nvalid,
            this->current_context->get_device(device));
      }
    }

  }
  // normalize images :
  // get the sum value per subap

  if (this->nstreams > 1) {
    subap_reduce_async(this->npix * this->npix, this->nvalid, this->streams,
        this->d_bincube->getData(), this->d_subsum->getData());
  } else {
    subap_reduce(this->d_bincube->getNbElem(), this->npix * this->npix,
        this->nvalid, this->d_bincube->getData(), this->d_subsum->getData(), this->current_context->get_device(device));
  }

  if (this->nstreams > 1) {
    subap_norm_async(this->d_bincube->getData(), this->d_bincube->getData(),
        this->d_fluxPerSub->getData(), this->d_subsum->getData(), this->nphot,
        this->npix * this->npix, this->d_bincube->getNbElem(), this->streams,
        current_context->get_device(device));
  } else {
    // multiply each subap by nphot*fluxPersub/sumPerSub
    subap_norm(this->d_bincube->getData(), this->d_bincube->getData(),
        this->d_fluxPerSub->getData(), this->d_subsum->getData(), this->nphot,
        this->npix * this->npix, this->d_bincube->getNbElem(), current_context->get_device(device));
  }
  //fprintf(stderr, "[%s@%d]: I'm here!\n", __FILE__, __LINE__);

  // add noise
  if (this->noise > -1) {
    //cout << "adding poisson noise" << endl;
	  this->d_bincube->prng('P');
  }
  if (this->noise > 0) {
    //cout << "adding detector noise" << endl;
    this->d_bincube->prng('N', this->noise, 1.0f);
  }
  return EXIT_SUCCESS;
}

//////////////////////////////
// PYRAMID WAVEFRONT SENSOR //
//////////////////////////////

// It starts by looking for the type of sensor. By default it assumes
// a pyramid wfs. The pyramid can also be explicitely asked for, or
// a roof prism can be asked for as well.

int sutra_wfs::comp_pyr_generic() {
    
    //START COMMENTING HERE TO SWITCH TO PYRAMID
  
  pyr_getpup(this->d_camplipup->getData(),
      this->d_gs->d_phase->d_screen->getData(), this->d_phalfxy->getData(),
      this->d_pupil->getData(), this->ntot, this->current_context->get_device(device));

  carma_fft(this->d_camplipup->getData(), this->d_camplifoc->getData(), -1,
      *this->d_camplipup->getPlan());
//
  pyr_submask(this->d_camplifoc->getData(), this->d_submask->getData(),
      this->ntot, this->current_context->get_device(device));

  cutilSafeCall(
      cudaMemset(this->d_hrimg->getData(), 0,
          sizeof(float) * this->d_hrimg->getNbElem()));
//
  //this->npup = 1;
  for (int cpt = 0; cpt < this->npup; cpt++) {
//    // modulation loop
//    // computes the high resolution images
    cutilSafeCall(
        cudaMemset(this->d_fttotim->getData(), 0,
            sizeof(cuFloatComplex) * this->d_fttotim->getNbElem()));
//
//    // here we split the image in 4 quadrant and roll them
    pyr_rollmod(this->d_fttotim->getData(), this->d_camplifoc->getData(),
        this->d_poffsets->getData(), (this->pyr_cx->getData())[cpt],
        (this->pyr_cy->getData())[cpt], this->ntot, this->nfft, this->current_context->get_device(device));
//
//    // case of diffractive pyramid
//    // multiply d_camplifoc->getData() by pyramid + modulation phase
//    // fft
//    // reorder the 4 quadrants
//
//    /*
//     pyr_rollmod(this->d_fttotim->getData(),this->d_camplifoc->getData(), this->d_poffsets->getData(),0,
//     0,this->ntot , this->nfft, this->current_context->get_device(device));
//     */

    carma_fft(this->d_fttotim->getData(), this->d_fttotim->getData(), 1,
        *this->d_fttotim->getPlan());

    float fact = 1.0f / this->nfft / this->nfft / this->nfft / 2.0;
    //if (cpt == this->npup-1) fact = fact / this->npup;

    pyr_abs2(this->d_hrimg->getData(), this->d_fttotim->getData(), fact,
        this->nfft, 4, this->current_context->get_device(device));
  }
//  /*
//   // spatial filtering by the pixel extent:
//   carma_fft(this->d_fttotim->getData(), this->d_fttotim->getData(), -1,
//   *this->d_fttotim->getPlan());
//
//   pyr_submask3d(this->d_fttotim->getData(), this->d_sincar->getData(),this->nfft, 4, this->current_context->get_device(device));
//
//   carma_fft(this->d_fttotim->getData(), this->d_fttotim->getData(), 1,
//   *this->d_fttotim->getPlan());
//
//   pyr_abs(this->d_hrimg->getData(), this->d_fttotim->getData(),this->nfft, 4, this->current_context->get_device(device));
//
//  pyr_fact(this->d_hrimg->getData(),1.0f/this->nfft/this->nfft,this->nfft,4,this->current_context->get_device(device));
//   */

    cutilSafeCall(
        cudaMemset(this->d_bincube->getData(), 0,
           sizeof(float) * this->d_bincube->getNbElem()));
           
 pyr_fillbin(this->d_bincube->getData(), this->d_hrimg->getData(),
      this->nrebin, this->nfft, this->nfft / this->nrebin, 4, this->current_context->get_device(device));

 pyr_subsum(this->d_subsum->getData(), this->d_bincube->getData(),
     this->d_validsubsx->getData(), this->d_validsubsy->getData(),
    this->nfft / this->nrebin, this->nvalid, 4, this->current_context->get_device(device));

 int blocks, threads;
 getNumBlocksAndThreads(current_context->get_device(device), this->nvalid, blocks, threads);
 reduce(this->nvalid, threads, blocks, this->d_subsum->getData(),
     this->d_subsum->getData());

  pyr_fact(this->d_bincube->getData(), this->nphot, this->d_subsum->getData(),
     this->nfft / this->nrebin, 4, this->current_context->get_device(device));

  // add noise
  if (this->noise > -1) {
    //cout << "adding poisson noise" << endl;
	  this->d_bincube->prng('P');
  }
  if (this->noise > 0) {
    //cout << "adding detector noise" << endl;
    this->d_bincube->prng('N', this->noise, 1.0f);
  }
  
  pyr_subsum(this->d_subsum->getData(), this->d_bincube->getData(),
     this->d_validsubsx->getData(), this->d_validsubsy->getData(),
     this->nfft / this->nrebin, this->nvalid, 4, this->current_context->get_device(device));
//  /*
//  reduce(this->nvalid, threads, blocks, this->d_subsum->getData(),
//      this->d_subsum->getData());
//  */
  return EXIT_SUCCESS;

  //___________________________________________________________________
  // MODIF ROOF SENSOR

  //PYR_GETPUP: reads pupil & phase and computes rolled electric field
/*   pyr_getpup(this->d_camplipup->getData(),
       this->d_gs->d_phase->d_screen->getData(), this->d_phalfxy->getData(),
       this->d_pupil->getData(), this->ntot, this->current_context->get_device(device));

   carma_fft(this->d_camplipup->getData(), this->d_camplifoc->getData(), -1,
       *this->d_camplipup->getPlan());

   pyr_submask(this->d_camplifoc->getData(), this->d_submask->getData(),
       this->ntot, this->current_context->get_device(device));

   cutilSafeCall(
       cudaMemset(this->d_hrimg->getData(), 0,
           sizeof(float) * this->d_hrimg->getNbElem()));

   //this->npup = 1;
   for (int cpt = 0; cpt < this->npup; cpt++) {
     // modulation loop
     // computes the high resolution images
     cutilSafeCall(
         cudaMemset(this->d_fttotim->getData(), 0,
             sizeof(cuFloatComplex) * this->d_fttotim->getNbElem()));

     roof_rollmod(this->d_fttotim->getData(), this->d_camplifoc->getData(),
         this->d_poffsets->getData(), (this->pyr_cx->getData())[cpt],
         (this->pyr_cy->getData())[cpt], this->ntot, this->nfft, this->current_context->get_device(device));
  //   
  //    pyr_rollmod(this->d_fttotim->getData(),this->d_camplifoc->getData(), this->d_poffsets->getData(),0,
  //    0,this->ntot , this->nfft, this->current_context->get_device(device));
  //    

     carma_fft(this->d_fttotim->getData(), this->d_fttotim->getData(), 1,
         *this->d_fttotim->getPlan());

     float fact = 1.0f / this->nfft / this->nfft / this->nfft / 2.0;
     //if (cpt == this->npup-1) fact = fact / this->npup;

     roof_abs2(this->d_hrimg->getData(), this->d_fttotim->getData(), fact,
         this->nfft, 4, this->current_context->get_device(device));
   }

   if (this->noise > 0) {
     this->d_bincube->prng('N', this->noise);
   } else
     cutilSafeCall(
         cudaMemset(this->d_bincube->getData(), 0,
             sizeof(float) * this->d_bincube->getNbElem()));

   roof_fillbin(this->d_bincube->getData(), this->d_hrimg->getData(),
       this->nrebin, this->nfft, this->nfft / this->nrebin, 4, this->current_context->get_device(device));

   pyr_subsum(this->d_subsum->getData(), this->d_bincube->getData(),
       this->d_validsubsx->getData(), this->d_validsubsy->getData(),
       this->nfft / this->nrebin, this->nvalid, 4, this->current_context->get_device(device));

   int blocks, threads;
   getNumBlocksAndThreads(this->current_context->get_device(device), this->nvalid, blocks, threads);
   reduce(this->nvalid, threads, blocks, this->d_subsum->getData(),
       this->d_subsum->getData());

   pyr_fact(this->d_bincube->getData(), this->nphot, this->d_subsum->getData(),
       this->nfft / this->nrebin, 4, this->current_context->get_device(device));

   pyr_subsum(this->d_subsum->getData(), this->d_bincube->getData(),
       this->d_validsubsx->getData(), this->d_validsubsy->getData(),
       this->nfft / this->nrebin, this->nvalid, 4, this->current_context->get_device(device));

   return EXIT_SUCCESS;
*/
}

int sutra_wfs::comp_image() {

  int result;
  if (this->type == "sh"){
    result = comp_sh_generic();
    if(result==EXIT_SUCCESS) {
      if (noise > 0) this->d_binimg->prng('N',this->noise);
      fillbinimg(this->d_binimg->getData(),this->d_bincube->getData(),this->npix,this->nvalid,this->npix*this->nxsub,
          this->d_validsubsx->getData(),this->d_validsubsy->getData(), 0/*(this->noise > 0)*/ ,this->current_context->get_device(device));
    }
  } else if (this->type == "pyr"){
    result = comp_pyr_generic();
  } else {
    DEBUG_TRACE("unknown wfs type : %s\n", this->type.c_str());
    result = EXIT_FAILURE;
  }
  return result;
}

int sutra_wfs::comp_image_tele() {
  int result;
  if (this->type == "sh"){
    result = comp_sh_generic();
  }else {
    DEBUG_TRACE("unknown wfs type : %s\n", this->type.c_str());
    result = EXIT_FAILURE;
  }

  if (noise > 0) this->d_binimg->prng('N',this->noise);

  if(result==EXIT_SUCCESS)
       fillbinimg_async(this->image_telemetry, this->d_binimg->getData(),
           this->d_bincube->getData(), this->npix, this->nvalid,
           this->npix * this->nxsub, this->d_validsubsx->getData(),
           this->d_validsubsy->getData(), this->d_binimg->getNbElem(), false,
           this->current_context->get_device(device));

  return result;
}

int sutra_wfs::slopes_geom(int type, float *slopes) {
  /*
   normalization notes :
   σ² = 0.17 (λ/D)^2 (D/r_0)^(5/3) , σ² en radians d'angle
   σ = sqrt(0.17 (λ/D)^2 (D/r_0)^(5/3)) * 206265 , σ en secondes

   // computing subaperture phase difference at edges

   todo : integrale( x * phase ) / integrale (x^2);
   with x = span(-0.5,0.5,npixels)(,-:1:npixels) * subap_diam * 2 * pi / lambda / 0.206265
   */
  if (type == 0) {
    // this is to convert in arcsec
    //> 206265* 0.000001/ 2 / 3.14159265 = 0.0328281
    // it would have been the case if the phase was given in radiants
    // but it is given in microns so normalization factor is
    // just 206265* 0.000001 = 0.206265

    //float alpha = 0.0328281 * this->d_gs->lambda / this->subapd;
    float alpha = 0.206265 / this->subapd;
    phase_reduce(this->nphase, this->nvalid,
        this->d_gs->d_phase->d_screen->getData(), slopes,
        this->d_phasemap->getData(), alpha);
  }

  if (type == 1) {

    //float alpha = 0.0328281 * this->d_gs->lambda / this->subapd;
    float alpha = 0.206265 / this->subapd;
    phase_derive(this->nphase * this->nphase * this->nvalid,
        this->nphase * this->nphase, this->nvalid, this->nphase,
        this->d_gs->d_phase->d_screen->getData(), slopes,
        this->d_phasemap->getData(), this->d_pupil->getData(), alpha,
        this->d_fluxPerSub->getData());

  }

  return EXIT_SUCCESS;
}

int sutra_wfs::slopes_geom(int type) {
  this->slopes_geom(type, this->d_slopes->getData());

  return EXIT_SUCCESS;
}

sutra_sensors::sutra_sensors(carma_context *context, const char* type, int nwfs,
    long *nxsub, long *nvalid, long *npix, long *nphase, long *nrebin,
    long *nfft, long *ntot, long npup, float *pdiam, float *nphot, int *lgs,
    int device) {

  for (int i = 0; i < nwfs; i++) {
    d_wfs.push_back(
        new sutra_wfs(context, type, nxsub[i], nvalid[i], npix[i], nphase[i],
            nrebin[i], nfft[i], ntot[i], npup, pdiam[i], nphot[i], lgs[i],
            device));
  }
}

sutra_sensors::sutra_sensors(carma_context *context, int nwfs, long *nxsub,
    long *nvalid, long *nphase, long npup, float *pdiam, int device) {

  for (int i = 0; i < nwfs; i++) {
    d_wfs.push_back(
        new sutra_wfs(context, nxsub[i], nvalid[i], nphase[i], npup, pdiam[i],
            device));
  }
}

sutra_sensors::~sutra_sensors() {
//  for (size_t idx = 0; idx < (this->d_wfs).size(); idx++) {
  while((this->d_wfs).size()>0) {
    delete this->d_wfs.back();
    d_wfs.pop_back();
  }
}

int sutra_sensors::sensors_initgs(float *xpos, float *ypos, float *lambda,
    float *mag, long *size, float *noise, long *seed) {
  for (size_t idx = 0; idx < (this->d_wfs).size(); idx++) {
    (this->d_wfs)[idx]->wfs_initgs(xpos[idx], ypos[idx], lambda[idx], mag[idx],
        size[idx], noise[idx], seed[idx]);
  }
  return EXIT_SUCCESS;
}
int sutra_sensors::sensors_initgs(float *xpos, float *ypos, float *lambda,
    float *mag, long *size, float *noise) {
  for (size_t idx = 0; idx < (this->d_wfs).size(); idx++) {
    (this->d_wfs)[idx]->wfs_initgs(xpos[idx], ypos[idx], lambda[idx], mag[idx],
        size[idx], noise[idx], 1234 * idx);
  }
  return EXIT_SUCCESS;
}
int sutra_sensors::sensors_initgs(float *xpos, float *ypos, float *lambda,
    float *mag, long *size) {
  for (size_t idx = 0; idx < (this->d_wfs).size(); idx++) {
    (this->d_wfs)[idx]->wfs_initgs(xpos[idx], ypos[idx], lambda[idx], mag[idx],
        size[idx], -1, 1234);
  }
  return EXIT_SUCCESS;
}

