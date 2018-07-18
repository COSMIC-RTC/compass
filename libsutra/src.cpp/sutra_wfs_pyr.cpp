#include <carma_utils.h>
#include <sutra_utils.h>
#include <sutra_wfs_pyr.h>

sutra_wfs_pyr::sutra_wfs_pyr(carma_context *context, sutra_telescope *d_tel,
                             carma_obj<cuFloatComplex> *d_camplipup,
                             carma_obj<cuFloatComplex> *d_camplifoc,
                             carma_obj<cuFloatComplex> *d_fttotim, long nxsub,
                             long nvalid, long npix, long nphase, long nrebin,
                             long nfft, long ntot, long npup, float pdiam,
                             float nphotons, float nphot4imat, int lgs,
                             bool roket, int device, const char *type_pyr)
    : sutra_wfs(context, d_tel, d_camplipup, d_camplifoc, d_fttotim, type_pyr,
                nxsub, nvalid, npix, nphase, nrebin, nfft, ntot, npup, pdiam,
                nphotons, nphot4imat, lgs, roket, device) {
  context->set_activeDevice(device, 1);

  long dims_data1[2];
  dims_data1[0] = 1;
  long dims_data2[3];
  dims_data2[0] = 2;
  long dims_data3[4];
  dims_data3[0] = 3;

  if (strcmp(type_pyr, "pyrhr") == 0) {
    dims_data2[1] = nfft;
    dims_data2[2] = nfft;

    this->d_hrimg =
        new carma_obj<float>(context, dims_data2);  // Useless for SH
    this->d_submask = new carma_obj<float>(context, dims_data2);
    this->d_camplipup = new carma_obj<cuFloatComplex>(context, dims_data2);
    this->d_camplifoc = new carma_obj<cuFloatComplex>(context, dims_data2);
    cufftHandle *plan = this->d_camplipup->getPlan();  ///< FFT plan
    carmafftSafeCall(
        cufftPlan2d(plan, dims_data2[1], dims_data2[2], CUFFT_C2C));

    this->d_fttotim = new carma_obj<cuFloatComplex>(context, dims_data2);

    this->d_poffsets = new carma_obj<cuFloatComplex>(context, dims_data2);
    this->d_phalfxy = new carma_obj<cuFloatComplex>(context, dims_data2);
    this->d_sincar = new carma_obj<float>(context, dims_data2);

    dims_data2[1] = nphase * nphase;
    dims_data2[2] = nvalid;
    this->d_phasemap = new carma_obj<int>(current_context, dims_data2);

    dims_data2[1] = nfft / nrebin;
    dims_data2[2] = nfft / nrebin;

    this->d_binimg = new carma_obj<float>(context, dims_data2);

    dims_data1[1] = 1;
    this->d_psum = new carma_obj<float>(context, dims_data1);

    if (this->roket) {
      this->d_binimg_notnoisy = new carma_obj<float>(context, dims_data2);
    }
    // using 1 stream for telemetry
    this->image_telemetry =
        new carma_host_obj<float>(dims_data2, MA_PAGELOCK, 1);
    this->nstreams = 1;
    while (nvalid % this->nstreams != 0) nstreams--;
    // err << "wfs uses " << nstreams << " streams" << endl;
    this->streams = new carma_streams(nstreams);

    dims_data1[1] = 2 * nvalid;
    this->d_slopes = new carma_obj<float>(context, dims_data1);

    dims_data2[1] = nfft;
    dims_data2[2] = nfft;
    dims_data3[3] = 4;

    dims_data1[1] = npup;
    this->pyr_cx = new carma_host_obj<float>(dims_data1, MA_WRICOMB);
    this->pyr_cy = new carma_host_obj<float>(dims_data1, MA_WRICOMB);

    dims_data1[1] = nvalid;
    this->d_subsum = new carma_obj<float>(context, dims_data1);

    this->d_fluxPerSub = new carma_obj<float>(context, dims_data1);
    dims_data1[1] = nvalid * 4;
    this->d_validsubsx = new carma_obj<int>(context, dims_data1);
    this->d_validsubsy = new carma_obj<int>(context, dims_data1);

  } else {  // pyr or roof
    dims_data3[1] = nfft;
    dims_data3[2] = nfft;
    dims_data3[3] = 4;

    this->d_hrimg =
        new carma_obj<float>(context, dims_data3);  // Useless for SH

    int mdims[2];

    dims_data2[1] = ntot;
    dims_data2[2] = ntot;
    this->d_submask = new carma_obj<float>(context, dims_data2);
    this->d_camplipup = new carma_obj<cuFloatComplex>(context, dims_data2);
    this->d_camplifoc = new carma_obj<cuFloatComplex>(context, dims_data2);
    cufftHandle *plan = this->d_camplipup->getPlan();  ///< FFT plan
    carmafftSafeCall(
        cufftPlan2d(plan, dims_data2[1], dims_data2[2], CUFFT_C2C));

    dims_data3[1] = nfft / nrebin;
    dims_data3[2] = nfft / nrebin;

    this->d_bincube = new carma_obj<float>(context, dims_data3);

    dims_data2[1] = nfft / nrebin * 2 + 3;
    dims_data2[2] = nfft / nrebin * 2 + 3;

    this->d_binimg = new carma_obj<float>(context, dims_data2);

    // using 1 stream for telemetry
    this->image_telemetry =
        new carma_host_obj<float>(dims_data2, MA_PAGELOCK, 1);

    this->nstreams = 1;
    while (nvalid % this->nstreams != 0) nstreams--;
    std::cerr << "wfs uses " << nstreams << " streams" << std::endl;
    this->streams = new carma_streams(nstreams);

    dims_data1[1] = 2 * nvalid;
    this->d_slopes = new carma_obj<float>(context, dims_data1);

    dims_data3[1] = nfft;
    dims_data3[2] = nfft;
    dims_data3[3] = 4;
    this->d_fttotim = new carma_obj<cuFloatComplex>(context, dims_data3);
    mdims[0] = (int)dims_data3[1];
    mdims[1] = (int)dims_data3[2];
    plan = this->d_fttotim->getPlan();  ///< FFT plan
    carmafftSafeCall(cufftPlanMany(plan, 2, mdims, NULL, 1, 0, NULL, 1, 0,
                                   CUFFT_C2C, (int)dims_data3[3]));

    //  dims_data2[1] = ntot;
    //  dims_data2[2] = ntot;
    //  this->d_pupil = new carma_obj<float>(context, dims_data2);

    dims_data1[1] = npup;
    this->pyr_cx = new carma_host_obj<float>(dims_data1, MA_WRICOMB);
    this->pyr_cy = new carma_host_obj<float>(dims_data1, MA_WRICOMB);
    dims_data2[1] = nfft;
    dims_data2[2] = nfft;
    this->d_poffsets = new carma_obj<cuFloatComplex>(context, dims_data2);
    dims_data2[1] = ntot;
    dims_data2[2] = ntot;
    this->d_phalfxy = new carma_obj<cuFloatComplex>(context, dims_data2);
    dims_data2[1] = nfft;
    dims_data2[2] = nfft;
    this->d_sincar = new carma_obj<float>(context, dims_data2);

    // dims_data1[1] = nvalid;
    // this->d_subsum = new carma_obj<float>(context, dims_data1);

    this->d_psum = new carma_obj<float>(context, dims_data1);

    this->d_fluxPerSub = new carma_obj<float>(context, dims_data1);
    this->d_validsubsx = new carma_obj<int>(context, dims_data1);
    this->d_validsubsy = new carma_obj<int>(context, dims_data1);

    dims_data2[1] = nphase * nphase;
    dims_data2[2] = nvalid;
    this->d_phasemap = new carma_obj<int>(context, dims_data2);
  }
}

sutra_wfs_pyr::~sutra_wfs_pyr() {
  current_context->set_activeDevice(device, 1);
  if (this->d_camplipup != 0L) delete this->d_camplipup;
  if (this->d_camplifoc != 0L) delete this->d_camplifoc;

  if (this->d_fttotim != 0L) delete this->d_fttotim;

  if (this->d_ftkernel != 0L) delete this->d_ftkernel;

  if (this->d_hrimg != 0L) delete this->d_hrimg;
  if (this->d_bincube != 0L) delete this->d_bincube;
  if (this->d_binimg != 0L) delete this->d_binimg;
  if (this->d_subsum != 0L) delete this->d_subsum;
  if (this->d_offsets != 0L) delete this->d_offsets;
  if (this->d_fluxPerSub != 0L) delete this->d_fluxPerSub;
  if (this->d_sincar != 0L) delete this->d_sincar;
  if (this->d_submask != 0L) delete this->d_submask;
  if (this->d_hrmap != 0L) delete this->d_hrmap;

  if (this->d_slopes != 0L) delete this->d_slopes;

  if (this->image_telemetry != 0L) delete this->image_telemetry;

  if (this->d_phasemap != 0L) delete this->d_phasemap;
  if (this->d_validsubsx != 0L) delete this->d_validsubsx;
  if (this->d_validsubsy != 0L) delete this->d_validsubsy;

  if (this->d_psum != 0L) delete this->d_psum;
  if (this->d_phalfxy != 0L) delete this->d_phalfxy;
  if (this->d_poffsets != 0L) delete this->d_poffsets;
  if (this->pyr_cx != 0L) delete this->pyr_cx;
  if (this->pyr_cy != 0L) delete this->pyr_cy;

  if (this->lgs) delete this->d_gs->d_lgs;

  if (this->d_gs != 0L) delete this->d_gs;

  delete this->streams;

  // delete this->current_context;
}

int sutra_wfs_pyr::loadarrays(cuFloatComplex *halfxy, cuFloatComplex *offsets,
                              float *focmask, float *cx, float *cy,
                              float *sincar, int *phasemap, int *validsubsx,
                              int *validsubsy) {
  current_context->set_activeDevice(device, 1);
  this->d_phalfxy->host2device(halfxy);
  this->d_poffsets->host2device(offsets);
  this->d_submask->host2device(focmask);
  this->pyr_cx->fill_from(cx);
  this->pyr_cy->fill_from(cy);
  this->d_sincar->host2device(sincar);
  this->d_phasemap->host2device(phasemap);
  this->d_validsubsx->host2device(validsubsx);
  this->d_validsubsy->host2device(validsubsy);

  return EXIT_SUCCESS;
}

int sutra_wfs_pyr::fill_binimage(int async) {
  if (this->d_binimg == NULL) {
    DEBUG_TRACE(
        "ERROR : d_bincube not initialized, did you do the allocate_buffers?");
    throw "ERROR : d_bincube not initialized, did you do the allocate_buffers?";
  }
  if (noise > 0) this->d_binimg->prng('N', this->noise);

  this->current_context->set_activeDevice(device, 1);
  if (async) {
    //    fillbinimg_async(this->image_telemetry, this->d_binimg->getData(),
    //        this->d_bincube->getData(), this->npix, this->nvalid_tot,
    //        this->npix * this->nxsub, this->d_validsubsx->getData(),
    //        this->d_validsubsy->getData(), this->d_binimg->getNbElem(), false,
    //        this->current_context->get_device(device));
    DEBUG_TRACE("ERROR : async version of fill_binimage not implemented...");
    throw "ERROR : async version of fill_binimage not implemented...";
  } else {
    pyr_fillbinimg(this->d_binimg->getData(), this->d_bincube->getData(),
                   this->nfft / this->nrebin, false,
                   this->current_context->get_device(device));
  }

  return EXIT_SUCCESS;
}
