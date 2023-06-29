// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_wfs_sh.cpp
//! \ingroup   libsutra
//! \class     SutraWfsSH
//! \brief     this class provides the wfs_sh features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.4
//! \date      2022/01/24

#include <carma_utils.h>
#include <sutra_utils.h>
#include <sutra_wfs_sh.h>
#include <cmath>

SutraWfsSH::SutraWfsSH(CarmaContext *context, SutraTelescope *d_tel,
                           CarmaObj<cuFloatComplex> *d_camplipup,
                           CarmaObj<cuFloatComplex> *d_camplifoc,
                           CarmaObj<cuFloatComplex> *d_fttotim, long nxsub,
                           long nvalid, long npix, long nphase, long nrebin,
                           long nfft, long ntot, long npup, float pdiam,
                           float nphotons, float nphot4imat, int lgs,
                           bool fakecam, int max_flux_per_pix, int max_pix_value,
                           bool is_low_order, bool roket, int device)
    : SutraWfs(context, d_tel, d_camplipup, d_camplifoc, d_fttotim, "sh",
                nxsub, nvalid, npix, nphase, nrebin, nfft, ntot, npup, pdiam,
                nphotons, nphot4imat, lgs, fakecam, max_flux_per_pix, max_pix_value,
                is_low_order, roket, device),
      d_binmap(nullptr),
      d_validpuppixx(nullptr),
      d_validpuppixy(nullptr),
      d_fsamplipup(nullptr),
      d_fsamplifoc(nullptr),
      fsampli_plan(nullptr) {}

int SutraWfsSH::define_mpi_rank(int rank, int size) {
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

int SutraWfsSH::allocate_buffers(
    map<vector<int>, cufftHandle *> campli_plans,
    map<vector<int>, cufftHandle *> fttotim_plans) {
  current_context->set_active_device(device, 1);
  long *dims_data1 = new long[2];
  dims_data1[0] = 1;
  long *dims_data2 = new long[3];
  dims_data2[0] = 2;
  long *dims_data3 = new long[4];
  dims_data3[0] = 3;

  dims_data2[1] = npix * nxsub;
  dims_data2[2] = npix * nxsub;

  if (rank == 0) {
    this->d_binimg = new CarmaObj<float>(current_context, dims_data2);
    if (this->roket) {
      this->d_binimg_notnoisy =
          new CarmaObj<float>(current_context, dims_data2);
    }
    if (this->fakecam) {
      this->d_camimg = new CarmaObj<uint16_t>(current_context, dims_data2);
      this->d_dark = new CarmaObj<float>(current_context, dims_data2);
      this->d_dark->reset();
      this->d_flat = new CarmaObj<float>(current_context, dims_data2);
      this->d_flat->memset(1);
    }
    // using 1 stream for telemetry
    this->image_telemetry =
        new CarmaHostObj<float>(dims_data2, MA_PAGELOCK, 1);
  }

  dims_data3[1] = nfft;
  dims_data3[2] = nfft;

  dims_data3[3] = nvalid;
  int mdims[2];

  // this->d_submask = new CarmaObj<float>(current_context, dims_data3); //
  // Useless for SH
  if (this->d_camplipup == nullptr)
    this->d_camplipup =
        new CarmaObj<cuFloatComplex>(current_context, dims_data3);
  if (this->d_camplifoc == nullptr)
    this->d_camplifoc =
        new CarmaObj<cuFloatComplex>(current_context, dims_data3);
  mdims[0] = (int)dims_data3[1];
  mdims[1] = (int)dims_data3[2];

  int vector_dims[3] = {mdims[0], mdims[1], (int)dims_data3[3]};
  vector<int> vdims(vector_dims,
                    vector_dims + sizeof(vector_dims) / sizeof(int));
  // vector<int> vdims(dims_data3 + 1, dims_data3 + 4);

  if (campli_plans.find(vdims) == campli_plans.end()) {
    // DEBUG_TRACE("Creating FFT plan : %d %d
    // %d",mdims[0],mdims[1],dims_data3[3]);print_mem_info();
    cufftHandle *plan = (cufftHandle *)malloc(
        sizeof(cufftHandle));  // = this->d_camplipup->get_plan(); ///< FFT plan
    carmafft_safe_call(cufftPlanMany(plan, 2, mdims, NULL, 1, 0, NULL, 1, 0,
                                   CUFFT_C2C, (int)dims_data3[3]));

    campli_plans.insert(pair<vector<int>, cufftHandle *>(vdims, plan));

    this->campli_plan = plan;
    // DEBUG_TRACE("FFT plan created");print_mem_info();
  } else {
    // DEBUG_TRACE("FFT plan already exists : %d %d
    // %d",mdims[0],mdims[1],dims_data3[3]);
    this->campli_plan = campli_plans.at(vdims);
  }

  dims_data3[1] = npix;
  dims_data3[2] = npix;
  if (rank == 0) dims_data3[3] = nvalid_tot;

  this->d_bincube = new CarmaObj<float>(current_context, dims_data3);

  this->nstreams = 1;
  while (nvalid % this->nstreams != 0) nstreams--;
  // std::cerr << "wfs uses " << nstreams << " streams" << std::endl;
  this->streams = new CarmaStreams(nstreams);

  dims_data1[1] = 2 * nvalid;
  this->d_slopes = new CarmaObj<float>(current_context, dims_data1);

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

    if (this->d_fttotim == nullptr)
      this->d_fttotim =
          new CarmaObj<cuFloatComplex>(current_context, dims_data3);

    mdims[0] = (int)dims_data3[1];
    mdims[1] = (int)dims_data3[2];
    int vector_dims[3] = {mdims[0], mdims[1], (int)dims_data3[3]};
    vector<int> vdims(vector_dims,
                      vector_dims + sizeof(vector_dims) / sizeof(int));

    if (fttotim_plans.find(vdims) == fttotim_plans.end()) {
      cufftHandle *plan = (cufftHandle *)malloc(
          sizeof(cufftHandle));  // = this->d_fttotim->get_plan(); ///< FFT plan
      // DEBUG_TRACE("Creating FFT plan :%d %d
      // %d",mdims[0],mdims[1],dims_data3[3]);print_mem_info();
      carmafft_safe_call(cufftPlanMany(plan, 2, mdims, NULL, 1, 0, NULL, 1, 0,
                                     CUFFT_C2C, (int)dims_data3[3]));
      fttotim_plans.insert(pair<vector<int>, cufftHandle *>(vdims, plan));
      this->fttotim_plan = plan;
      // DEBUG_TRACE("FFT plan created : ");print_mem_info();
    } else {
      // DEBUG_TRACE("FFT plan already exists %d %d
      // %d",mdims[0],mdims[1],dims_data3[3]);
      this->fttotim_plan = fttotim_plans.at(vdims);
    }
    dims_data1[1] = nfft * nfft;
    this->d_hrmap = new CarmaObj<int>(current_context, dims_data1);

  } else {
    if (this->lgs) {
      dims_data3[1] = ntot;
      dims_data3[2] = ntot;
      dims_data3[3] = nvalid;
      // this->d_fttotim = new CarmaObj<cuFloatComplex>(current_context,
      // dims_data3);
      mdims[0] = (int)dims_data3[1];
      mdims[1] = (int)dims_data3[2];
      int vector_dims[3] = {mdims[0], mdims[1], (int)dims_data3[3]};
      vector<int> vdims(vector_dims,
                        vector_dims + sizeof(vector_dims) / sizeof(int));

      if (fttotim_plans.find(vdims) == fttotim_plans.end()) {
        // DEBUG_TRACE("Creating FFT plan : %d %d
        // %d",mdims[0],mdims[1],dims_data3[3]);print_mem_info();
        cufftHandle *plan = (cufftHandle *)malloc(sizeof(
            cufftHandle));  // = this->d_fttotim->get_plan(); ///< FFT plan
        carmafft_safe_call(cufftPlanMany(plan, 2, mdims, NULL, 1, 0, NULL, 1, 0,
                                       CUFFT_C2C, (int)dims_data3[3]));
        fttotim_plans.insert(pair<vector<int>, cufftHandle *>(vdims, plan));
        this->fttotim_plan = plan;
        // DEBUG_TRACE("FFT plan created : ");print_mem_info();
      } else {
        // DEBUG_TRACE("FFT plan already exists : %d %d
        // %d",mdims[0],mdims[1],dims_data3[3]);
        this->fttotim_plan = fttotim_plans.at(vdims);
      }
    }
  }

  dims_data2[1] = ntot;
  dims_data2[2] = ntot;
  this->d_ftkernel = new CarmaObj<cuFloatComplex>(current_context, dims_data2);

  dims_data2[1] = nphase;
  dims_data2[2] = nphase;
  this->d_offsets = new CarmaObj<float>(current_context, dims_data2);

  dims_data2[1] = nrebin * nrebin;
  dims_data2[2] = npix * npix;
  this->d_binmap = new CarmaObj<int>(current_context, dims_data2);

  dims_data1[1] = nvalid_tot;
  this->d_intensities = new CarmaObj<float>(current_context, dims_data1);

  this->d_fluxPerSub = new CarmaObj<float>(current_context, dims_data1);
  this->d_validsubsx = new CarmaObj<int>(current_context, dims_data1);
  this->d_validsubsy = new CarmaObj<int>(current_context, dims_data1);
  this->d_validpuppixx = new CarmaObj<int>(current_context, dims_data1);
  this->d_validpuppixy = new CarmaObj<int>(current_context, dims_data1);

  dims_data2[1] = nphase * nphase;
  dims_data2[2] = nvalid;
  this->d_phasemap = new CarmaObj<int>(current_context, dims_data2);
  
  dims_data2[1] = nphase * nphase;
  dims_data2[2] = nvalid * 2;
  this->d_ttprojmat = new CarmaObj<float>(current_context, dims_data2);

  dims_data1[1] = nphase * nphase * nvalid;
  this->d_ttprojvec = new CarmaObj<float>(current_context, dims_data1);

  delete[] dims_data1;
  delete[] dims_data2;
  delete[] dims_data3;
  return EXIT_SUCCESS;
}

SutraWfsSH::~SutraWfsSH() {
  current_context->set_active_device(device, 1);

  if (this->d_ftkernel != 0L) delete this->d_ftkernel;

  if (this->d_bincube != 0L) delete this->d_bincube;
  if (this->d_binimg != 0L) delete this->d_binimg;
  if (this->d_intensities != 0L) delete this->d_intensities;
  if (this->d_offsets != 0L) delete this->d_offsets;
  if (this->d_fluxPerSub != 0L) delete this->d_fluxPerSub;
  if (this->d_sincar != 0L) delete this->d_sincar;
  if (this->d_hrmap != 0L) delete this->d_hrmap;

  if (this->d_slopes != 0L) delete this->d_slopes;

  if (this->image_telemetry != 0L) delete this->image_telemetry;

  if (this->d_phasemap != 0L) delete this->d_phasemap;
  if (this->d_ttprojmat != 0L) delete this->d_ttprojmat;
  if (this->d_ttprojvec != 0L) delete this->d_ttprojvec;
  if (this->d_binmap != 0L) delete this->d_binmap;
  if (this->d_validsubsx != 0L) delete this->d_validsubsx;
  if (this->d_validsubsy != 0L) delete this->d_validsubsy;
  if (this->d_validpuppixx != 0L) delete this->d_validpuppixx;
  if (this->d_validpuppixy != 0L) delete this->d_validpuppixy;

  if (this->lgs) delete this->d_gs->d_lgs;

  if (this->d_gs != 0L) delete this->d_gs;

  delete this->streams;

  // delete this->current_context;
}

int SutraWfsSH::load_arrays(int *phasemap, int *hrmap, int *binmap,
                             float *offsets, float *fluxPerSub, int *validsubsx,
                             int *validsubsy, int *istart, int *jstart, 
                             float *ttprojmat, cuFloatComplex *kernel) {
  if (this->d_bincube == NULL) {
    DEBUG_TRACE(
        "ERROR : d_bincube not initialized, did you do the allocate_buffers?");
    throw "ERROR : d_bincube not initialized, did you do the allocate_buffers?";
  }
  current_context->set_active_device(device, 1);
  this->d_phasemap->host2device(&phasemap[offset * nphase * nphase]);
  this->d_ttprojmat->host2device(ttprojmat);
  this->d_offsets->host2device(offsets);
  this->d_binmap->host2device(binmap);
  this->d_fluxPerSub->host2device(&fluxPerSub[offset]);
  if (this->ntot != this->nfft) this->d_hrmap->host2device(hrmap);
  this->d_validsubsx->host2device(&validsubsx[offset]);
  this->d_validsubsy->host2device(&validsubsy[offset]);
  this->d_validpuppixx->host2device(istart);
  this->d_validpuppixy->host2device(jstart);
  this->d_ftkernel->host2device(kernel);

  return EXIT_SUCCESS;
}

/////////////////////////////////////////////////////////
// COMPUTATION OF THE SHACK-HARTMANN WAVEFRONT SENSOR  //
/////////////////////////////////////////////////////////

int SutraWfsSH::comp_generic() {
  if (this->d_bincube == NULL) {
    DEBUG_TRACE(
        "ERROR : d_bincube not initialized, did you do the allocate_buffers?");
    throw "ERROR : d_bincube not initialized, did you do the allocate_buffers?";
  }
  current_context->set_active_device(device, 1);

  carma_safe_call(
      cudaMemset(this->d_camplipup->get_data(), 0,
                 sizeof(cuFloatComplex) * this->d_camplipup->get_nb_elements()));
  /*
   carma_safe_call(
   cudaMemset(this->d_camplifoc->get_data(), 0,
   sizeof(cuFloatComplex) * this->d_camplifoc->get_nb_elements()));
   */

  // segment phase and fill cube of complex ampli with exp(i*phase_seg)
  if (this->d_fsamplipup == nullptr) {// No Field stop case
    fillcamplipup(
        this->d_camplipup->get_data(), this->d_gs->d_phase->d_screen->get_data(),
        this->d_offsets->get_data(), this->d_pupil->get_data(), this->d_gs->scale,
        this->d_validpuppixx->get_data(), this->d_validpuppixy->get_data(),
        this->d_validsubsx->get_data(), this->d_validsubsy->get_data(),
        this->nphase, this->d_gs->d_phase->d_screen->get_dims(1), this->nfft,
        this->nphase * this->nphase * this->nvalid,
        this->current_context->get_device(device), 0);  //, this->offset);
                                                        // do fft of the cube
  }
  else { // There is a field stop
    cudaMemset(this->d_fsamplipup->get_data(), 0,
                sizeof(cuFloatComplex) * this->d_fsamplipup->get_nb_elements());
    fillfsamplipup(
        this->d_fsamplipup->get_data(), this->d_gs->d_phase->d_screen->get_data(),
        this->d_pupil->get_data(), this->d_gs->scale, 
        this->d_fsamplipup->get_dims(1),this->d_gs->d_phase->d_screen->get_dims(1), 
        this->current_context->get_device(device));

    CarmaFFT(this->d_fsamplipup->get_data(), this->d_fsamplifoc->get_data(), 1,
            *this->fsampli_plan);
    apply_submask(this->d_fsamplifoc->get_data(),
               this->d_submask->get_data(),
               this->d_fsamplipup->get_dims(1), this->current_context->get_device(device));
    CarmaFFT(this->d_fsamplifoc->get_data(), this->d_fsamplifoc->get_data(), -1,
            *this->fsampli_plan);
    fillcamplipup(
        this->d_camplipup->get_data(), this->d_fsamplifoc->get_data(),
        this->d_offsets->get_data(),
        this->d_validpuppixx->get_data(), this->d_validpuppixy->get_data(),
        this->d_validsubsx->get_data(), this->d_validsubsy->get_data(),
        this->nphase, this->d_fsamplifoc->get_dims(1), this->nfft,
        this->nphase * this->nphase * this->nvalid,
        this->current_context->get_device(device));  //, this->offset);

  }

  CarmaFFT(this->d_camplipup->get_data(), this->d_camplifoc->get_data(), 1,
            *this->campli_plan);  //*this->d_camplipup->get_plan());

  // get the hrimage by taking the | |^2
  // keep it in amplifoc to save mem space
  abs2c(this->d_camplifoc->get_data(), this->d_camplifoc->get_data(),
        this->nfft * this->nfft * this->nvalid,
        current_context->get_device(device));

  // set bincube to 0 or noise
  carma_safe_call(cudaMemset(this->d_bincube->get_data(), 0,
                           sizeof(float) * this->d_bincube->get_nb_elements()));
  // increase fov if required
  // and fill bincube with data from hrimg
  // we need to do this sequentially if nvalid > nmaxhr to
  // keep raesonable mem occupancy
  if (this->ntot != this->nfft) {
    for (int cc = 0; cc < this->nffthr; cc++) {
      carma_safe_call(
          cudaMemset(this->d_fttotim->get_data(), 0,
                     sizeof(cuFloatComplex) * this->d_fttotim->get_nb_elements()));

      int indxstart1, indxstart2 = 0, indxstart3;

      if ((cc == this->nffthr - 1) && (this->nvalid % this->nmaxhr != 0)) {
        indxstart1 = (this->nfft * this->nfft * this->nvalid) -
                     this->nfft * this->nfft * this->nmaxhr;
        if (this->lgs)
          indxstart2 = this->ntot * this->nvalid - this->ntot * this->nmaxhr;
        indxstart3 = this->d_bincube->get_nb_elements() -
                     this->npix * this->npix * this->nmaxhr;
      } else {
        indxstart1 = this->nfft * this->nfft * this->nmaxhr * cc;
        if (this->lgs) indxstart2 = this->ntot * this->nmaxhr * cc;
        indxstart3 = this->npix * this->npix * this->nmaxhr * cc;
      }

      cuFloatComplex *data = this->d_camplifoc->get_data();
      indexfill(this->d_fttotim->get_data(), &(data[indxstart1]),
                this->d_hrmap->get_data(), this->nfft, this->ntot,
                this->nfft * this->nfft * this->nmaxhr,
                this->current_context->get_device(device));

      if (this->lgs) {
        // compute lgs spot on the fly from binned profile image
        this->d_gs->d_lgs->lgs_makespot(
            this->current_context->get_device(device), indxstart2);
        // convolve with psf
        CarmaFFT(this->d_fttotim->get_data(), this->d_fttotim->get_data(), 1,
                  *this->fttotim_plan);  //*this->d_fttotim->get_plan());

        convolve(this->d_fttotim->get_data(),
                 this->d_gs->d_lgs->d_ftlgskern->get_data(),
                 this->ntot * this->ntot * this->nmaxhr,
                 this->current_context->get_device(device));

        CarmaFFT(this->d_fttotim->get_data(), this->d_fttotim->get_data(), -1,
                  *this->fttotim_plan);  //*this->d_fttotim->get_plan());
      }

      if (this->kernconv) {
        CarmaFFT(this->d_fttotim->get_data(), this->d_fttotim->get_data(), 1,
                  *this->fttotim_plan);  //*this->d_fttotim->get_plan());

        convolve_cube(this->d_fttotim->get_data(), this->d_ftkernel->get_data(),
                      this->ntot * this->ntot * this->nmaxhr,
                      this->d_ftkernel->get_nb_elements(),
                      this->current_context->get_device(device));

        CarmaFFT(this->d_fttotim->get_data(), this->d_fttotim->get_data(), -1,
                  *this->fttotim_plan);  //*this->d_fttotim->get_plan());
      }

      float *data2 = this->d_bincube->get_data();
      // fprintf(stderr, "[%s@%d]: I'm here!\n", __FILE__, __LINE__);
      if (this->nstreams > 1) {
        fillbincube_async(this->streams, &(data2[indxstart3]),
                          this->d_fttotim->get_data(), this->d_binmap->get_data(),
                          this->ntot * this->ntot, this->npix * this->npix,
                          this->nrebin * this->nrebin, this->nmaxhr,
                          this->current_context->get_device(device));
      } else {
        fillbincube(&(data2[indxstart3]), this->d_fttotim->get_data(),
                    this->d_binmap->get_data(), this->ntot * this->ntot,
                    this->npix * this->npix, this->nrebin * this->nrebin,
                    this->nmaxhr, this->current_context->get_device(device));
        // fprintf(stderr, "[%s@%d]: I'm here!\n", __FILE__, __LINE__);
      }
    }
  } else {
    if (this->lgs) {
      carma_safe_call(
          cudaMemset(this->d_fttotim->get_data(), 0,
                     sizeof(cuFloatComplex) * this->d_fttotim->get_nb_elements()));
      this->d_gs->d_lgs->lgs_makespot(this->current_context->get_device(device),
                                      0);

      CarmaFFT(this->d_camplifoc->get_data(), this->d_fttotim->get_data(), 1,
                *this->fttotim_plan);  //*this->d_fttotim->get_plan());

      convolve(this->d_fttotim->get_data(),
               this->d_gs->d_lgs->d_ftlgskern->get_data(),
               this->ntot * this->ntot * this->nvalid,
               this->current_context->get_device(device));

      CarmaFFT(this->d_fttotim->get_data(), this->d_fttotim->get_data(), -1,
                *this->fttotim_plan);  //*this->d_fttotim->get_plan());

      if (this->nstreams > 1) {
        fillbincube_async(this->streams, this->d_bincube->get_data(),
                          this->d_fttotim->get_data(), this->d_binmap->get_data(),
                          this->nfft * this->nfft, this->npix * this->npix,
                          this->nrebin * this->nrebin, this->nvalid,
                          this->current_context->get_device(device));
      } else {
        fillbincube(this->d_bincube->get_data(), this->d_fttotim->get_data(),
                    this->d_binmap->get_data(), this->nfft * this->nfft,
                    this->npix * this->npix, this->nrebin * this->nrebin,
                    this->nvalid, this->current_context->get_device(device));
      }
    } else {
      if (this->kernconv) {
        CarmaFFT(this->d_camplifoc->get_data(), this->d_camplifoc->get_data(), 1,
                  *this->campli_plan);  //*this->d_camplipup->get_plan());

        convolve_cube(this->d_camplifoc->get_data(), this->d_ftkernel->get_data(),
                      this->nfft * this->nfft * this->nvalid,
                      this->d_ftkernel->get_nb_elements(),
                      this->current_context->get_device(device));

        CarmaFFT(this->d_camplifoc->get_data(), this->d_camplifoc->get_data(),
                  -1, *this->campli_plan);  //*this->d_camplipup->get_plan());
      }

      if (this->nstreams > 1) {
        fillbincube_async(this->streams, this->d_bincube->get_data(),
                          this->d_camplifoc->get_data(),
                          this->d_binmap->get_data(), this->nfft * this->nfft,
                          this->npix * this->npix, this->nrebin * this->nrebin,
                          this->nvalid,
                          this->current_context->get_device(device));
      } else {
        fillbincube(this->d_bincube->get_data(), this->d_camplifoc->get_data(),
                    this->d_binmap->get_data(), this->nfft * this->nfft,
                    this->npix * this->npix, this->nrebin * this->nrebin,
                    this->nvalid, this->current_context->get_device(device));
      }
    }
  }
  // normalize images :
  // get the sum value per subap

  if (this->nstreams > 1) {
    subap_reduce_async(this->npix * this->npix, this->nvalid, this->streams,
                       this->d_bincube->get_data(),
                       this->d_intensities->get_data());
  } else {
    subap_reduce(this->d_bincube->get_nb_elements(), this->npix * this->npix,
                 this->nvalid, this->d_bincube->get_data(),
                 this->d_intensities->get_data(),
                 this->current_context->get_device(device));
  }

  if (this->nstreams > 1) {
    subap_norm_async(this->d_bincube->get_data(), this->d_bincube->get_data(),
                     this->d_fluxPerSub->get_data(),
                     this->d_intensities->get_data(), this->nphot,
                     this->npix * this->npix, this->d_bincube->get_nb_elements(),
                     this->streams, current_context->get_device(device));
  } else {
    // multiply each subap by nphot*fluxPersub/sumPerSub
    subap_norm(this->d_bincube->get_data(), this->d_bincube->get_data(),
               this->d_fluxPerSub->get_data(), this->d_intensities->get_data(),
               this->nphot, this->npix * this->npix,
               this->d_bincube->get_nb_elements(),
               current_context->get_device(device));
  }
  // fprintf(stderr, "[%s@%d]: I'm here!\n", __FILE__, __LINE__);

  if (this->roket) {  // Get here the bincube before adding noise,
                      // usefull for error budget
    fillbinimg(this->d_binimg_notnoisy->get_data(), this->d_bincube->get_data(),
               this->npix, this->nvalid_tot, this->npix * this->nxsub,
               this->d_validsubsx->get_data(), this->d_validsubsy->get_data(),
               false, this->current_context->get_device(device));
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

int SutraWfsSH::fill_binimage(int async = 0) {
  if (this->d_binimg == NULL) {
    DEBUG_TRACE(
        "ERROR : d_bincube not initialized, did you do the allocate_buffers?");
    throw "ERROR : d_bincube not initialized, did you do the allocate_buffers?";
  }
  if (noise > 0) this->d_binimg->prng('N', this->noise);

  current_context->set_active_device(device, 1);
  if (async) {
    fillbinimg_async(this->image_telemetry, this->d_binimg->get_data(),
                     this->d_bincube->get_data(), this->npix, this->nvalid_tot,
                     this->npix * this->nxsub, this->d_validsubsx->get_data(),
                     this->d_validsubsy->get_data(), this->d_binimg->get_nb_elements(),
                     false, this->current_context->get_device(device));
  } else {
    fillbinimg(this->d_binimg->get_data(), this->d_bincube->get_data(),
               this->npix, this->nvalid_tot, this->npix * this->nxsub,
               this->d_validsubsx->get_data(), this->d_validsubsy->get_data(),
               false, this->current_context->get_device(device));
  }
  return EXIT_SUCCESS;
}

int SutraWfsSH::comp_image(bool noise) {
  current_context->set_active_device(device, 1);
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
  if (this->fakecam)
    result *= digitalize(this->d_camimg->get_data(), this->d_binimg->get_data(),
                         this->d_dark->get_data(), this->d_flat->get_data(),
                         this->max_flux_per_pix, this->max_pix_value,
                         this->d_binimg->get_nb_elements(),
                         this->current_context->get_device(this->device));

  return result;
}

int SutraWfsSH::comp_nphot(float ittime, float optthroughput, float diam,
                             int nxsub, float zerop, float gsmag,
                             float lgsreturnperwatt, float laserpower) {
  this->d_gs->mag = gsmag;
  if (laserpower == 0) {
    if (zerop == 0) zerop = 1e11;
    this->nphot = zerop * pow(10., (-0.4 * gsmag)) * ittime * optthroughput *
                  (diam / nxsub) * (diam / nxsub);
    // include throughput to WFS for unobstructed
    // subaperture per iteration
  } else {  // we are dealing with a LGS
    this->nphot = lgsreturnperwatt * laserpower * optthroughput *
                  (diam / nxsub) * (diam / nxsub) * 1e4 * ittime;
    // detected by WFS
    // ... for given power include throughput to WFS
    // for unobstructed subaperture per iteration
  }
  return EXIT_SUCCESS;
}

int SutraWfsSH::set_bincube(float *bincube, int nElem) {
  current_context->set_active_device(device, 1);
  if (nElem == this->d_bincube->get_nb_elements())
    this->d_bincube->host2device(bincube);
  else
    DEBUG_TRACE("Wrong size of cube");
  return EXIT_SUCCESS;
}

int SutraWfsSH::set_field_stop(map<vector<int>, cufftHandle *> campli_plans,
                                float* field_stop, int N) {
  if(this->d_submask != nullptr) {
    delete d_submask;
    delete d_fsamplipup;
    delete d_fsamplifoc;
    fsampli_plan = nullptr;
  }
  long dims_data[3] = {2, N, N};

  this->d_submask = new CarmaObj<float>(current_context, dims_data, field_stop);
  this->d_fsamplipup = new CarmaObj<cuFloatComplex>(current_context, dims_data);
  this->d_fsamplifoc = new CarmaObj<cuFloatComplex>(current_context, dims_data);

  vector<int> vector_dims {2, N, N};
  if (campli_plans.find(vector_dims) == campli_plans.end()) {
    // DEBUG_TRACE("Creating FFT plan : %d %d
    // %d",mdims[0],mdims[1],dims_data3[3]);print_mem_info();
    cufftHandle *plan = (cufftHandle *)malloc(
        sizeof(cufftHandle));  // = this->d_camplipup->get_plan(); ///< FFT plan
    carmafft_safe_call(cufftPlan2d(plan, N, N, CUFFT_C2C));

    campli_plans.insert(pair<vector<int>, cufftHandle *>(vector_dims, plan));

    this->fsampli_plan = plan;
    // DEBUG_TRACE("FFT plan created");print_mem_info();
  } else {
    // DEBUG_TRACE("FFT plan already exists : %d %d
    // %d",mdims[0],mdims[1],dims_data3[3]);
    this->fsampli_plan = campli_plans.at(vector_dims);
  }

  return EXIT_SUCCESS;
}
