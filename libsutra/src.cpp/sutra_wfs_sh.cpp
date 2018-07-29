#include <carma_utils.h>
#include <sutra_utils.h>
#include <sutra_wfs_sh.h>

sutra_wfs_sh::sutra_wfs_sh(carma_context *context, sutra_telescope *d_tel,
                           carma_obj<cuFloatComplex> *d_camplipup,
                           carma_obj<cuFloatComplex> *d_camplifoc,
                           carma_obj<cuFloatComplex> *d_fttotim, long nxsub,
                           long nvalid, long npix, long nphase, long nrebin,
                           long nfft, long ntot, long npup, float pdiam,
                           float nphotons, float nphot4imat, int lgs,
                           bool roket, int device)
    : sutra_wfs(context, d_tel, d_camplipup, d_camplifoc, d_fttotim, "sh",
                nxsub, nvalid, npix, nphase, nrebin, nfft, ntot, npup, pdiam,
                nphotons, nphot4imat, lgs, roket, device),
      d_binmap(nullptr),
      d_istart(nullptr),
      d_jstart(nullptr) {}

int sutra_wfs_sh::define_mpi_rank(int rank, int size) {
  if (this->device == 0) this->device = rank % current_context->get_ndevice();

  int count = this->nvalid / size;
  int i;
  this->rank = rank;
  int r = this->nvalid % size;
  if (rank < r) {
    this->nvalid = this->nvalid / size + 1;
    this->offset = rank * this->nvalid;
  } else {
    this->nvalid = this->nvalid / size;
    this->offset = rank * this->nvalid + r;
  }

  displ_bincube = new int[size];
  count_bincube = new int[size];
  for (i = 0; i < r; i++) {
    count_bincube[i] = npix * npix * (count + 1);
    displ_bincube[i] = i * count_bincube[i];
  }
  for (i = r; i < size; i++) {
    count_bincube[i] = npix * npix * count;
    displ_bincube[i] = i * count_bincube[i] + r * npix * npix;
  }

  return EXIT_SUCCESS;
}

int sutra_wfs_sh::allocate_buffers(
    map<vector<int>, cufftHandle *> campli_plans,
    map<vector<int>, cufftHandle *> fttotim_plans) {
  current_context->set_activeDevice(device, 1);
  long *dims_data1 = new long[2];
  dims_data1[0] = 1;
  long *dims_data2 = new long[3];
  dims_data2[0] = 2;
  long *dims_data3 = new long[4];
  dims_data3[0] = 3;

  dims_data2[1] = npix * nxsub;
  dims_data2[2] = npix * nxsub;

  if (rank == 0) {
    this->d_binimg = new carma_obj<float>(current_context, dims_data2);
    // using 1 stream for telemetry
    this->image_telemetry =
        new carma_host_obj<float>(dims_data2, MA_PAGELOCK, 1);
  }

  dims_data3[1] = nfft;
  dims_data3[2] = nfft;

  dims_data3[3] = nvalid;

  int mdims[2];

  // this->d_submask = new carma_obj<float>(current_context, dims_data3); //
  // Useless for SH

  // this->d_camplipup = new carma_obj<cuFloatComplex>(current_context,
  // dims_data3); this->d_camplifoc = new
  // carma_obj<cuFloatComplex>(current_context, dims_data3);
  mdims[0] = (int)dims_data3[1];
  mdims[1] = (int)dims_data3[2];

  int vector_dims[3] = {mdims[0], mdims[1], (int)dims_data3[3]};
  vector<int> vdims(vector_dims,
                    vector_dims + sizeof(vector_dims) / sizeof(int));
  // vector<int> vdims(dims_data3 + 1, dims_data3 + 4);

  if (campli_plans.find(vdims) == campli_plans.end()) {
    // DEBUG_TRACE("Creating FFT plan : %d %d
    // %d",mdims[0],mdims[1],dims_data3[3]);printMemInfo();
    cufftHandle *plan = (cufftHandle *)malloc(
        sizeof(cufftHandle));  // = this->d_camplipup->getPlan(); ///< FFT plan
    carmafftSafeCall(cufftPlanMany(plan, 2, mdims, NULL, 1, 0, NULL, 1, 0,
                                   CUFFT_C2C, (int)dims_data3[3]));

    campli_plans.insert(pair<vector<int>, cufftHandle *>(vdims, plan));

    this->campli_plan = plan;
    // DEBUG_TRACE("FFT plan created");printMemInfo();
  } else {
    // DEBUG_TRACE("FFT plan already exists : %d %d
    // %d",mdims[0],mdims[1],dims_data3[3]);
    this->campli_plan = campli_plans.at(vdims);
  }

  dims_data3[1] = npix;
  dims_data3[2] = npix;
  if (rank == 0) dims_data3[3] = nvalid_tot;

  this->d_bincube = new carma_obj<float>(current_context, dims_data3);
  if (this->roket) {
    this->d_bincube_notnoisy =
        new carma_obj<float>(current_context, dims_data3);
  }

  this->nstreams = 1;
  while (nvalid % this->nstreams != 0) nstreams--;
  std::cerr << "wfs uses " << nstreams << " streams" << std::endl;
  this->streams = new carma_streams(nstreams);

  dims_data1[1] = 2 * nvalid;
  this->d_slopes = new carma_obj<float>(current_context, dims_data1);

  if (this->ntot != this->nfft) {
    // this is the big array => we use nmaxhr and treat it sequentially
    int mnmax = 500;
    if (nvalid > 2 * mnmax) {
      nmaxhr = compute_nmaxhr(nvalid);

      this->nffthr = (nvalid % this->nmaxhr == 0 ? nvalid / this->nmaxhr
                                                 : nvalid / this->nmaxhr + 1);
    }
    dims_data3[1] = ntot;
    dims_data3[2] = ntot;
    dims_data3[3] = nmaxhr;

    // this->d_fttotim = new carma_obj<cuFloatComplex>(current_context,
    // dims_data3);

    mdims[0] = (int)dims_data3[1];
    mdims[1] = (int)dims_data3[2];
    int vector_dims[3] = {mdims[0], mdims[1], (int)dims_data3[3]};
    vector<int> vdims(vector_dims,
                      vector_dims + sizeof(vector_dims) / sizeof(int));

    if (fttotim_plans.find(vdims) == fttotim_plans.end()) {
      cufftHandle *plan = (cufftHandle *)malloc(
          sizeof(cufftHandle));  // = this->d_fttotim->getPlan(); ///< FFT plan
      // DEBUG_TRACE("Creating FFT plan :%d %d
      // %d",mdims[0],mdims[1],dims_data3[3]);printMemInfo();
      carmafftSafeCall(cufftPlanMany(plan, 2, mdims, NULL, 1, 0, NULL, 1, 0,
                                     CUFFT_C2C, (int)dims_data3[3]));
      fttotim_plans.insert(pair<vector<int>, cufftHandle *>(vdims, plan));
      this->fttotim_plan = plan;
      // DEBUG_TRACE("FFT plan created : ");printMemInfo();
    } else {
      // DEBUG_TRACE("FFT plan already exists %d %d
      // %d",mdims[0],mdims[1],dims_data3[3]);
      this->fttotim_plan = fttotim_plans.at(vdims);
    }
    dims_data1[1] = nfft * nfft;
    this->d_hrmap = new carma_obj<int>(current_context, dims_data1);

  } else {
    if (this->lgs) {
      dims_data3[1] = ntot;
      dims_data3[2] = ntot;
      dims_data3[3] = nvalid;
      // this->d_fttotim = new carma_obj<cuFloatComplex>(current_context,
      // dims_data3);
      mdims[0] = (int)dims_data3[1];
      mdims[1] = (int)dims_data3[2];
      int vector_dims[3] = {mdims[0], mdims[1], (int)dims_data3[3]};
      vector<int> vdims(vector_dims,
                        vector_dims + sizeof(vector_dims) / sizeof(int));

      if (fttotim_plans.find(vdims) == fttotim_plans.end()) {
        // DEBUG_TRACE("Creating FFT plan : %d %d
        // %d",mdims[0],mdims[1],dims_data3[3]);printMemInfo();
        cufftHandle *plan = (cufftHandle *)malloc(sizeof(
            cufftHandle));  // = this->d_fttotim->getPlan(); ///< FFT plan
        carmafftSafeCall(cufftPlanMany(plan, 2, mdims, NULL, 1, 0, NULL, 1, 0,
                                       CUFFT_C2C, (int)dims_data3[3]));
        fttotim_plans.insert(pair<vector<int>, cufftHandle *>(vdims, plan));
        this->fttotim_plan = plan;
        // DEBUG_TRACE("FFT plan created : ");printMemInfo();
      } else {
        // DEBUG_TRACE("FFT plan already exists : %d %d
        // %d",mdims[0],mdims[1],dims_data3[3]);
        this->fttotim_plan = fttotim_plans.at(vdims);
      }
    }
  }

  dims_data2[1] = ntot;
  dims_data2[2] = ntot;
  this->d_ftkernel = new carma_obj<cuFloatComplex>(current_context, dims_data2);

  dims_data2[1] = nphase;
  dims_data2[2] = nphase;
  this->d_offsets = new carma_obj<float>(current_context, dims_data2);

  dims_data1[1] = nxsub;
  this->d_istart = new carma_obj<int>(current_context, dims_data1);
  this->d_jstart = new carma_obj<int>(current_context, dims_data1);

  dims_data2[1] = nrebin * nrebin;
  dims_data2[2] = npix * npix;
  this->d_binmap = new carma_obj<int>(current_context, dims_data2);

  dims_data1[1] = nvalid_tot;
  this->d_subsum = new carma_obj<float>(current_context, dims_data1);

  this->d_fluxPerSub = new carma_obj<float>(current_context, dims_data1);
  this->d_validsubsx = new carma_obj<int>(current_context, dims_data1);
  this->d_validsubsy = new carma_obj<int>(current_context, dims_data1);

  dims_data2[1] = nphase * nphase;
  dims_data2[2] = nvalid;
  this->d_phasemap = new carma_obj<int>(current_context, dims_data2);

  delete[] dims_data1;
  delete[] dims_data2;
  delete[] dims_data3;
  return EXIT_SUCCESS;
}

sutra_wfs_sh::~sutra_wfs_sh() {
  current_context->set_activeDevice(device, 1);

  if (this->d_ftkernel != 0L) delete this->d_ftkernel;

  if (this->d_bincube != 0L) delete this->d_bincube;
  if (this->d_binimg != 0L) delete this->d_binimg;
  if (this->d_subsum != 0L) delete this->d_subsum;
  if (this->d_offsets != 0L) delete this->d_offsets;
  if (this->d_fluxPerSub != 0L) delete this->d_fluxPerSub;
  if (this->d_sincar != 0L) delete this->d_sincar;
  if (this->d_hrmap != 0L) delete this->d_hrmap;

  if (this->d_slopes != 0L) delete this->d_slopes;

  if (this->image_telemetry != 0L) delete this->image_telemetry;

  if (this->d_phasemap != 0L) delete this->d_phasemap;
  if (this->d_binmap != 0L) delete this->d_binmap;
  if (this->d_validsubsx != 0L) delete this->d_validsubsx;
  if (this->d_validsubsy != 0L) delete this->d_validsubsy;
  if (this->d_istart != 0L) delete this->d_istart;
  if (this->d_jstart != 0L) delete this->d_jstart;

  if (this->lgs) delete this->d_gs->d_lgs;

  if (this->d_gs != 0L) delete this->d_gs;

  delete this->streams;

  // delete this->current_context;
}

int sutra_wfs_sh::loadarrays(int *phasemap, int *hrmap, int *binmap,
                             float *offsets, float *fluxPerSub, int *validsubsx,
                             int *validsubsy, int *istart, int *jstart,
                             cuFloatComplex *kernel) {
  if (this->d_bincube == NULL) {
    DEBUG_TRACE(
        "ERROR : d_bincube not initialized, did you do the allocate_buffers?");
    throw "ERROR : d_bincube not initialized, did you do the allocate_buffers?";
  }
  current_context->set_activeDevice(device, 1);
  this->d_phasemap->host2device(&phasemap[offset * nphase * nphase]);
  this->d_offsets->host2device(offsets);
  this->d_binmap->host2device(binmap);
  this->d_fluxPerSub->host2device(&fluxPerSub[offset]);
  if (this->ntot != this->nfft) this->d_hrmap->host2device(hrmap);
  this->d_validsubsx->host2device(&validsubsx[offset]);
  this->d_validsubsy->host2device(&validsubsy[offset]);
  this->d_istart->host2device(istart);
  this->d_jstart->host2device(jstart);
  this->d_ftkernel->host2device(kernel);

  return EXIT_SUCCESS;
}

/////////////////////////////////////////////////////////
// COMPUTATION OF THE SHACK-HARTMANN WAVEFRONT SENSOR  //
/////////////////////////////////////////////////////////

int sutra_wfs_sh::comp_generic() {
  if (this->d_bincube == NULL) {
    DEBUG_TRACE(
        "ERROR : d_bincube not initialized, did you do the allocate_buffers?");
    throw "ERROR : d_bincube not initialized, did you do the allocate_buffers?";
  }
  current_context->set_activeDevice(device, 1);

  carmaSafeCall(
      cudaMemset(this->d_camplipup->getData(), 0,
                 sizeof(cuFloatComplex) * this->d_camplipup->getNbElem()));
  /*
   carmaSafeCall(
   cudaMemset(this->d_camplifoc->getData(), 0,
   sizeof(cuFloatComplex) * this->d_camplifoc->getNbElem()));
   */

  // segment phase and fill cube of complex ampli with exp(i*phase_seg)
  fillcamplipup(
      this->d_camplipup->getData(), this->d_gs->d_phase->d_screen->getData(),
      this->d_offsets->getData(), this->d_pupil->getData(), this->d_gs->scale,
      this->d_istart->getData(), this->d_jstart->getData(),
      this->d_validsubsx->getData(), this->d_validsubsy->getData(),
      this->nphase, this->d_gs->d_phase->d_screen->getDims(1), this->nfft,
      this->nphase * this->nphase * this->nvalid,
      this->current_context->get_device(device), 0);  //, this->offset);
                                                      // do fft of the cube
  carma_fft(this->d_camplipup->getData(), this->d_camplifoc->getData(), 1,
            *this->campli_plan);  //*this->d_camplipup->getPlan());

  // get the hrimage by taking the | |^2
  // keep it in amplifoc to save mem space
  abs2c(this->d_camplifoc->getData(), this->d_camplifoc->getData(),
        this->nfft * this->nfft * this->nvalid,
        current_context->get_device(device));

  // set bincube to 0 or noise
  carmaSafeCall(cudaMemset(this->d_bincube->getData(), 0,
                           sizeof(float) * this->d_bincube->getNbElem()));
  // increase fov if required
  // and fill bincube with data from hrimg
  // we need to do this sequentially if nvalid > nmaxhr to
  // keep raesonable mem occupancy
  if (this->ntot != this->nfft) {
    for (int cc = 0; cc < this->nffthr; cc++) {
      carmaSafeCall(
          cudaMemset(this->d_fttotim->getData(), 0,
                     sizeof(cuFloatComplex) * this->d_fttotim->getNbElem()));

      int indxstart1, indxstart2 = 0, indxstart3;

      if ((cc == this->nffthr - 1) && (this->nvalid % this->nmaxhr != 0)) {
        indxstart1 = (this->nfft * this->nfft * this->nvalid) -
                     this->nfft * this->nfft * this->nmaxhr;
        if (this->lgs)
          indxstart2 = this->ntot * this->nvalid - this->ntot * this->nmaxhr;
        indxstart3 = this->d_bincube->getNbElem() -
                     this->npix * this->npix * this->nmaxhr;
      } else {
        indxstart1 = this->nfft * this->nfft * this->nmaxhr * cc;
        if (this->lgs) indxstart2 = this->ntot * this->nmaxhr * cc;
        indxstart3 = this->npix * this->npix * this->nmaxhr * cc;
      }

      cuFloatComplex *data = this->d_camplifoc->getData();
      indexfill(this->d_fttotim->getData(), &(data[indxstart1]),
                this->d_hrmap->getData(), this->nfft, this->ntot,
                this->nfft * this->nfft * this->nmaxhr,
                this->current_context->get_device(device));

      if (this->lgs) {
        // compute lgs spot on the fly from binned profile image
        this->d_gs->d_lgs->lgs_makespot(
            this->current_context->get_device(device), indxstart2);
        // convolve with psf
        carma_fft(this->d_fttotim->getData(), this->d_fttotim->getData(), 1,
                  *this->fttotim_plan);  //*this->d_fttotim->getPlan());

        convolve(this->d_fttotim->getData(),
                 this->d_gs->d_lgs->d_ftlgskern->getData(),
                 this->ntot * this->ntot * this->nmaxhr,
                 this->current_context->get_device(device));

        carma_fft(this->d_fttotim->getData(), this->d_fttotim->getData(), -1,
                  *this->fttotim_plan);  //*this->d_fttotim->getPlan());
      }

      if (this->kernconv) {
        carma_fft(this->d_fttotim->getData(), this->d_fttotim->getData(), 1,
                  *this->fttotim_plan);  //*this->d_fttotim->getPlan());

        convolve_cube(this->d_fttotim->getData(), this->d_ftkernel->getData(),
                      this->ntot * this->ntot * this->nmaxhr,
                      this->d_ftkernel->getNbElem(),
                      this->current_context->get_device(device));

        carma_fft(this->d_fttotim->getData(), this->d_fttotim->getData(), -1,
                  *this->fttotim_plan);  //*this->d_fttotim->getPlan());
      }

      float *data2 = this->d_bincube->getData();
      // fprintf(stderr, "[%s@%d]: I'm here!\n", __FILE__, __LINE__);
      if (this->nstreams > 1) {
        fillbincube_async(this->streams, &(data2[indxstart3]),
                          this->d_fttotim->getData(), this->d_binmap->getData(),
                          this->ntot * this->ntot, this->npix * this->npix,
                          this->nrebin * this->nrebin, this->nmaxhr,
                          this->current_context->get_device(device));
      } else {
        fillbincube(&(data2[indxstart3]), this->d_fttotim->getData(),
                    this->d_binmap->getData(), this->ntot * this->ntot,
                    this->npix * this->npix, this->nrebin * this->nrebin,
                    this->nmaxhr, this->current_context->get_device(device));
        // fprintf(stderr, "[%s@%d]: I'm here!\n", __FILE__, __LINE__);
      }
    }
  } else {
    if (this->lgs) {
      carmaSafeCall(
          cudaMemset(this->d_fttotim->getData(), 0,
                     sizeof(cuFloatComplex) * this->d_fttotim->getNbElem()));
      this->d_gs->d_lgs->lgs_makespot(this->current_context->get_device(device),
                                      0);

      carma_fft(this->d_camplifoc->getData(), this->d_fttotim->getData(), 1,
                *this->fttotim_plan);  //*this->d_fttotim->getPlan());

      convolve(this->d_fttotim->getData(),
               this->d_gs->d_lgs->d_ftlgskern->getData(),
               this->ntot * this->ntot * this->nvalid,
               this->current_context->get_device(device));

      carma_fft(this->d_fttotim->getData(), this->d_fttotim->getData(), -1,
                *this->fttotim_plan);  //*this->d_fttotim->getPlan());

      if (this->nstreams > 1) {
        fillbincube_async(this->streams, this->d_bincube->getData(),
                          this->d_fttotim->getData(), this->d_binmap->getData(),
                          this->nfft * this->nfft, this->npix * this->npix,
                          this->nrebin * this->nrebin, this->nvalid,
                          this->current_context->get_device(device));
      } else {
        fillbincube(this->d_bincube->getData(), this->d_fttotim->getData(),
                    this->d_binmap->getData(), this->nfft * this->nfft,
                    this->npix * this->npix, this->nrebin * this->nrebin,
                    this->nvalid, this->current_context->get_device(device));
      }
    } else {
      if (this->kernconv) {
        carma_fft(this->d_camplifoc->getData(), this->d_camplifoc->getData(), 1,
                  *this->campli_plan);  //*this->d_camplipup->getPlan());

        convolve_cube(this->d_camplifoc->getData(), this->d_ftkernel->getData(),
                      this->nfft * this->nfft * this->nvalid,
                      this->d_ftkernel->getNbElem(),
                      this->current_context->get_device(device));

        carma_fft(this->d_camplifoc->getData(), this->d_camplifoc->getData(),
                  -1, *this->campli_plan);  //*this->d_camplipup->getPlan());
      }

      if (this->nstreams > 1) {
        fillbincube_async(this->streams, this->d_bincube->getData(),
                          this->d_camplifoc->getData(),
                          this->d_binmap->getData(), this->nfft * this->nfft,
                          this->npix * this->npix, this->nrebin * this->nrebin,
                          this->nvalid,
                          this->current_context->get_device(device));
      } else {
        fillbincube(this->d_bincube->getData(), this->d_camplifoc->getData(),
                    this->d_binmap->getData(), this->nfft * this->nfft,
                    this->npix * this->npix, this->nrebin * this->nrebin,
                    this->nvalid, this->current_context->get_device(device));
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
                 this->nvalid, this->d_bincube->getData(),
                 this->d_subsum->getData(),
                 this->current_context->get_device(device));
  }

  if (this->nstreams > 1) {
    subap_norm_async(this->d_bincube->getData(), this->d_bincube->getData(),
                     this->d_fluxPerSub->getData(), this->d_subsum->getData(),
                     this->nphot, this->npix * this->npix,
                     this->d_bincube->getNbElem(), this->streams,
                     current_context->get_device(device));
  } else {
    // multiply each subap by nphot*fluxPersub/sumPerSub
    subap_norm(this->d_bincube->getData(), this->d_bincube->getData(),
               this->d_fluxPerSub->getData(), this->d_subsum->getData(),
               this->nphot, this->npix * this->npix,
               this->d_bincube->getNbElem(),
               current_context->get_device(device));
  }
  // fprintf(stderr, "[%s@%d]: I'm here!\n", __FILE__, __LINE__);

  if (this->roket) {  // Get here the bincube before adding noise,
                      // usefull for error budget
    this->d_bincube->copyInto(this->d_bincube_notnoisy->getData(),
                              this->d_bincube->getNbElem());
  }
  // add noise
  if (this->noise > -1) {
    // cout << "adding poisson noise" << endl;
    this->d_bincube->prng('P');
  }
  if (this->noise > 0) {
    // cout << "adding detector noise" << endl;
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

int sutra_wfs_sh::fill_binimage(int async = 0) {
  if (this->d_binimg == NULL) {
    DEBUG_TRACE(
        "ERROR : d_bincube not initialized, did you do the allocate_buffers?");
    throw "ERROR : d_bincube not initialized, did you do the allocate_buffers?";
  }
  if (noise > 0) this->d_binimg->prng('N', this->noise);

  current_context->set_activeDevice(device, 1);
  if (async) {
    fillbinimg_async(this->image_telemetry, this->d_binimg->getData(),
                     this->d_bincube->getData(), this->npix, this->nvalid_tot,
                     this->npix * this->nxsub, this->d_validsubsx->getData(),
                     this->d_validsubsy->getData(), this->d_binimg->getNbElem(),
                     false, this->current_context->get_device(device));
  } else {
    fillbinimg(this->d_binimg->getData(), this->d_bincube->getData(),
               this->npix, this->nvalid_tot, this->npix * this->nxsub,
               this->d_validsubsx->getData(), this->d_validsubsy->getData(),
               false, this->current_context->get_device(device));
  }
  return EXIT_SUCCESS;
}

int sutra_wfs_sh::comp_image(bool noise) {
  current_context->set_activeDevice(device, 1);
  int result;
  if (noise)
    result = comp_generic();
  else {
    float tmp = this->noise;
    this->noise = -1.0;
    result = comp_generic();
    this->noise = tmp;
  }
  result *= this->fill_binimage();

  return result;
}
