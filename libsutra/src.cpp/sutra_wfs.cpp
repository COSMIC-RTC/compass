// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_wfs.cpp
//! \ingroup   libsutra
//! \class     SutraWfs
//! \brief     this class provides the wfs features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.0
//! \date      2022/01/24

#include <sutra_wfs.h>

SutraWfs::SutraWfs(CarmaContext *context, SutraTelescope *d_tel,
                     CarmaObj<cuFloatComplex> *d_camplipup,
                     CarmaObj<cuFloatComplex> *d_camplifoc,
                     CarmaObj<cuFloatComplex> *d_fttotim, string type,
                     long nxsub, long nvalid, long npix, long nphase,
                     long nrebin, long nfft, long ntot, long npup, float pdiam,
                     float nphotons, float nphot4imat, int lgs, bool fakecam,
                     int max_flux_per_pix, int max_pix_value, bool is_low_order,
                     bool roket, int device)
    : device(device),
      type(type),
      nxsub(nxsub),
      nvalid(nvalid),
      npix(npix),
      nrebin(nrebin),
      nfft(nfft),
      ntot(ntot),
      npup(npup),
      nphase(nphase),
      nmaxhr(nvalid),
      nffthr(1),
      subapd(pdiam),
      nphot(nphotons),
      nphot4imat(nphot4imat),
      noise(0),
      lgs(lgs),
      fakecam(fakecam),
      max_flux_per_pix(max_flux_per_pix),
      max_pix_value(max_pix_value),
      is_low_order(is_low_order),
      kernconv(false),
      roket(roket),
      d_dark(nullptr),
      d_flat(nullptr),
      d_camimg(nullptr),
      campli_plan(nullptr),
      fttotim_plan(nullptr),
      d_ftkernel(nullptr),
      d_camplipup(nullptr),
      d_camplifoc(nullptr),
      d_fttotim(nullptr),
      d_pupil(d_tel->d_pupil_m),
      d_bincube(nullptr),
      d_binimg(nullptr),
      d_binimg_notnoisy(nullptr),
      d_intensities(nullptr),
      d_offsets(nullptr),
      d_fluxPerSub(nullptr),
      d_sincar(nullptr),
      d_hrmap(nullptr),
      d_slopes(nullptr),
      image_telemetry(nullptr),
      d_gs(nullptr),
      streams(nullptr),
      nstreams(0),
      d_phasemap(nullptr),
      d_validsubsx(nullptr),
      d_validsubsy(nullptr),
      d_ttprojmat(nullptr),
      d_ttprojvec(nullptr),
      current_context(context),
      offset(0),
      nvalid_tot(nvalid),
      rank(0),
      displ_bincube(nullptr),
      count_bincube(nullptr),
      d_submask(nullptr) {
  if (d_camplipup != nullptr && !is_low_order) this->d_camplipup = d_camplipup;
  if (d_camplifoc != nullptr && !is_low_order) this->d_camplifoc = d_camplifoc;
  if (d_fttotim != nullptr && !is_low_order) this->d_fttotim = d_fttotim;
}

int SutraWfs::wfs_initgs(CarmaObj<float> *d_lgskern,
                          CarmaObj<cuFloatComplex> *d_ftlgskern,
                          map<vector<int>, cufftHandle *> ftlgskern_plans,
                          float xpos, float ypos, float lambda, float mag,
                          float zerop, long size, float noise, long seed,
                          float G, float thetaML, float dx, float dy) {
  current_context->set_active_device(device, 1);
  this->d_gs = new SutraSource(current_context, xpos, ypos, lambda, mag, zerop,
                                size, "wfs", this->device);
  this->d_gs->G = G;
  this->d_gs->thetaML = thetaML;
  this->d_gs->dx = dx;
  this->d_gs->dy = dy;

  set_noise(noise, seed);

  if (this->lgs) {
    this->d_gs->d_lgs =
        new SutraLGS(current_context, d_lgskern, d_ftlgskern, ftlgskern_plans,
                      this->nvalid, this->ntot, this->nmaxhr);
    this->d_gs->lgs = this->lgs;
  }

  return EXIT_SUCCESS;
}

int SutraWfs::set_noise(float noise, long seed) {
  current_context->set_active_device(device, 1);
  this->noise = noise;
  if (this->type != "pyrhr") {
    if (noise > -1) {
      this->d_bincube->init_prng(seed);
      // this->d_bincube->prng('N', noise, 0.0f);
    }
    if (noise > 0) {
      this->d_binimg->init_prng(seed);
      // this->d_binimg->prng('N', noise, 0.0f);
    }
  } else {
    if (noise > -1) {
      this->d_binimg->init_prng(seed);
    }
  }
  return EXIT_SUCCESS;
}

int SutraWfs::set_pupil(float *pupil) {
  current_context->set_active_device(device, 1);
  this->d_pupil->host2device(pupil);
  int nbdevices = d_pupil_ngpu.size();
  for (int ndevice = 1; ndevice < nbdevices; ndevice++) {
    current_context->set_active_device(ndevice, 1);
    d_pupil_ngpu[ndevice]->host2device(pupil);
  }
  current_context->set_active_device(device, 1);
  return EXIT_SUCCESS;
}

int SutraWfs::load_kernels(float *lgskern) {
  current_context->set_active_device(device, 1);
  if (this->lgs)
    this->d_gs->d_lgs->load_kernels(lgskern,
                                    this->current_context->get_device(device));

  return EXIT_SUCCESS;
}

int SutraWfs::sensor_trace(SutraAtmos *yatmos) {
  current_context->set_active_device(device, 1);
  // do raytracing to get the phase
  this->d_gs->raytrace(yatmos);

  // same with textures
  // apparently not working for large sizes...
  // this->d_gs->raytrace_shm(yatmos);

  return EXIT_SUCCESS;
}

int SutraWfs::sensor_trace(SutraDms *ydm, int rst) {
  current_context->set_active_device(device, 1);
  // do raytracing to get the phase
  this->d_gs->raytrace(ydm, rst, 0);

  return EXIT_SUCCESS;
}

int SutraWfs::sensor_trace(int rst) {
  current_context->set_active_device(device, 1);
  // do raytracing to get the phase
  this->d_gs->raytrace(rst);

  return EXIT_SUCCESS;
}

int SutraWfs::sensor_trace(SutraAtmos *yatmos, SutraDms *ydms) {
  this->d_gs->raytrace(yatmos);
  this->d_gs->raytrace(ydms, 0, 0);

  return EXIT_SUCCESS;
}

int SutraWfs::slopes_geom(float *slopes, int type) {
  current_context->set_active_device(device, 1);
  /*
   normalization notes :
   ���� = 0.17 (��/D)^2 (D/r_0)^(5/3) , ���� en radians d'angle
   �� = sqrt(0.17 (��/D)^2 (D/r_0)^(5/3)) * 206265 , �� en secondes

   // computing subaperture phase difference at edges

   todo : integrale( x * phase ) / integrale (x^2);
   with x = span(-0.5,0.5,npixels)(,-:1:npixels) * subap_diam * 2 * pi / lambda
   / 0.206265
   */
  if (type == 0) {
    // this is to convert in arcsec
    //> 206265* 0.000001/ 2 / CARMA_PI = 0.0328281
    // it would have been the case if the phase was given in radiants
    // but it is given in microns so normalization factor is
    // just 206265* 0.000001 = 0.206265

    // float alpha = 0.0328281 * this->d_gs->lambda / this->subapd;
    float alpha = 0.206265 / this->subapd * this->nphase / (this->nphase - 1);
    phase_reduce(this->nphase, this->nvalid,
                 this->d_gs->d_phase->d_screen->get_data(), slopes,
                 this->d_phasemap->get_data(), alpha);
  }

  if (type == 1) {
    // float alpha = 0.0328281 * this->d_gs->lambda / this->subapd;
    float alpha = 0.206265 / this->subapd * this->nphase / (this->nphase - 1);
    phase_derive(this->nphase * this->nphase * this->nvalid,
                 this->nphase * this->nphase, this->nvalid, this->nphase,
                 this->d_gs->d_phase->d_screen->get_data(), slopes,
                 this->d_phasemap->get_data(), this->d_pupil->get_data(), alpha,
                 this->d_fluxPerSub->get_data());
  }
  
  if (type == 2) { // linear projection method
    phase_project(this->nphase, this->nvalid, 
                  this->d_gs->d_phase->d_screen->get_data(), slopes, 
                  this->d_phasemap->get_data(), this->d_ttprojmat->get_data(), 
                  this->d_ttprojvec->get_data(), 
                  current_context->get_device(device));
  }

  return EXIT_SUCCESS;
}

int SutraWfs::slopes_geom(int type) {
  this->slopes_geom(this->d_slopes->get_data(), type);

  return EXIT_SUCCESS;
}

int SutraWfs::set_binimg(float *binimg, int nElem) {
  current_context->set_active_device(device, 1);

  if (nElem == this->d_binimg->get_nb_elements()) {
    this->d_binimg->host2device(binimg);
  } else
    DEBUG_TRACE("Wrong dimension of binimg");
  return EXIT_SUCCESS;
}

int SutraWfs::set_dark(float *dark, int nElem) {
  current_context->set_active_device(device, 1);

  if (this->d_dark == nullptr) {
    std::cout << "WFS d_dark array has not been initialized yet" << std::endl;
    return EXIT_FAILURE;
  }

  if (nElem == this->d_dark->get_nb_elements()) {
    this->d_dark->host2device(dark);
  } else
    DEBUG_TRACE("Wrong dimension of dark");
  return EXIT_SUCCESS;
}

int SutraWfs::set_flat(float *flat, int nElem) {
  current_context->set_active_device(device, 1);

  if (this->d_flat == nullptr) {
    std::cout << "WFS d_flat array has not been initialized yet" << std::endl;
    return EXIT_FAILURE;
  }

  if (nElem == this->d_flat->get_nb_elements()) {
    this->d_flat->host2device(flat);
  } else
    DEBUG_TRACE("Wrong dimension of flat");
  return EXIT_SUCCESS;
}

int SutraWfs::set_fakecam(bool fakecam) {
  current_context->set_active_device(device, 1);

  this->fakecam = fakecam;
  if (this->d_camimg == nullptr) {
    this->d_camimg =
        new CarmaObj<uint16_t>(current_context, this->d_binimg->get_dims());
    this->d_dark =
        new CarmaObj<float>(current_context, this->d_binimg->get_dims());
    this->d_dark->reset();
    this->d_flat =
        new CarmaObj<float>(current_context, this->d_binimg->get_dims());
    this->d_flat->memset(1);
  }

  return EXIT_SUCCESS;
}

int SutraWfs::set_max_flux_per_pix(int max_flux_per_pix) {
  this->max_flux_per_pix = max_flux_per_pix;
  return EXIT_SUCCESS;
}

int SutraWfs::set_max_pix_value(int max_pix_value) {
  this->max_pix_value = max_pix_value;
  return EXIT_SUCCESS;
}
