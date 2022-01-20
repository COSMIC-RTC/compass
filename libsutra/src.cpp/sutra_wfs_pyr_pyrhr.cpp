// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2019 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.
//  Distributed under GNU - LGPL
//
//  COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
//  General Public License as published by the Free Software Foundation, either version 3 of the License,
//  or any later version.
//
//  COMPASS: End-to-end AO simulation tool using GPU acceleration
//  The COMPASS platform was designed to meet the need of high-performance for the simulation of AO systems.
//
//  The final product includes a software package for simulating all the critical subcomponents of AO,
//  particularly in the context of the ELT and a real-time core based on several control approaches,
//  with performances consistent with its integration into an instrument. Taking advantage of the specific
//  hardware architecture of the GPU, the COMPASS tool allows to achieve adequate execution speeds to
//  conduct large simulation campaigns called to the ELT.
//
//  The COMPASS platform can be used to carry a wide variety of simulations to both testspecific components
//  of AO of the E-ELT (such as wavefront analysis device with a pyramid or elongated Laser star), and
//  various systems configurations such as multi-conjugate AO.
//
//  COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
//  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//  See the GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License along with COMPASS.
//  If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>.
// -----------------------------------------------------------------------------

//! \file      sutra_wfs_pyr_pyrhr.cpp
//! \ingroup   libsutra
//! \class     SutraWfs_PyrHR
//! \brief     this class provides the wfs_pyr_pyrhr features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.2.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#include <carma_utils.h>
#include <sutra_utils.h>
#include <sutra_wfs_pyr_pyrhr.h>
#include <cmath>

SutraWfs_PyrHR::SutraWfs_PyrHR(
    CarmaContext *context, SutraTelescope *d_tel,
    CarmaObj<cuFloatComplex> *d_camplipup,
    CarmaObj<cuFloatComplex> *d_camplifoc,
    CarmaObj<cuFloatComplex> *d_fttotim, long nxsub, long nvalid, long npupils,
    long npix, long nphase, long nrebin, long nfft, long ntot, long npup,
    float pdiam, float nphotons, float nphot4imat, int lgs, bool fakecam,
    int max_flux_per_pix, int max_pix_value, bool roket, int device)
    : SutraWfs(context, d_tel, d_camplipup, d_camplifoc, d_fttotim, "pyrhr",
                nxsub, nvalid, npix, nphase, nrebin, nfft, ntot, npup, pdiam,
                nphotons, nphot4imat, lgs, fakecam, max_flux_per_pix, max_pix_value,
                false, roket, device) {
  context->set_active_device(device, 1);
  long dims_data1[2];
  dims_data1[0] = 1;
  long dims_data2[3];
  dims_data2[0] = 2;

  dims_data2[1] = nfft;
  dims_data2[2] = nfft;

  this->npupils = npupils;
  this->compute_pyrfocalplane = false;
  this->d_hrimg = new CarmaObj<float>(context, dims_data2);  // Useless for SH
  this->d_submask = new CarmaObj<float>(context, dims_data2);
  this->d_camplipup = new CarmaObj<cuFloatComplex>(context, dims_data2);
  this->d_camplifoc = new CarmaObj<cuFloatComplex>(context, dims_data2);
  this->d_pyrfocalplane = new CarmaObj<float>(context, dims_data2);
  cufftHandle *plan = this->d_camplipup->get_plan();  ///< FFT plan
  carmafft_safe_call(cufftPlan2d(plan, dims_data2[1], dims_data2[2], CUFFT_C2C));

  this->d_fttotim = new CarmaObj<cuFloatComplex>(context, dims_data2);

  this->d_poffsets = new CarmaObj<cuFloatComplex>(context, dims_data2);
  this->d_phalfxy = new CarmaObj<cuFloatComplex>(context, dims_data2);
  this->d_sincar = new CarmaObj<float>(context, dims_data2);

  dims_data2[1] = nphase * nphase;
  dims_data2[2] = nvalid;
  this->d_phasemap = new CarmaObj<int>(current_context, dims_data2);

  dims_data2[1] = nphase * nphase;
  dims_data2[2] = nvalid * 2;
  this->d_ttprojmat = new CarmaObj<float>(current_context, dims_data2);

  dims_data1[1] = nphase * nphase * nvalid;
  this->d_ttprojvec = new CarmaObj<float>(current_context, dims_data1);

  dims_data2[1] = nfft / nrebin;
  dims_data2[2] = nfft / nrebin;

  this->d_binimg = new CarmaObj<float>(context, dims_data2);
  if (this->fakecam) {
    this->d_camimg = new CarmaObj<uint16_t>(current_context, dims_data2);
    this->d_dark = new CarmaObj<float>(current_context, dims_data2);
    this->d_dark->reset();
    this->d_flat = new CarmaObj<float>(current_context, dims_data2);
    this->d_flat->memset(1);
  }

  dims_data1[1] = 1;
  this->d_psum = new CarmaObj<float>(context, dims_data1);

  if (this->roket) {
    this->d_binimg_notnoisy = new CarmaObj<float>(context, dims_data2);
  }
  // using 1 stream for telemetry
  this->image_telemetry = new CarmaHostObj<float>(dims_data2, MA_PAGELOCK, 1);
  this->nstreams = 1;
  while (nvalid % this->nstreams != 0) nstreams--;
  // err << "wfs uses " << nstreams << " streams" << endl;
  this->streams = new CarmaStreams(nstreams);

  dims_data1[1] = 2 * nvalid;
  this->d_slopes = new CarmaObj<float>(context, dims_data1);

  dims_data1[1] = npup;
  this->pyr_cx = new CarmaHostObj<float>(dims_data1, MA_WRICOMB);
  this->pyr_cy = new CarmaHostObj<float>(dims_data1, MA_WRICOMB);
  this->pyr_mod_weights = new CarmaHostObj<float>(dims_data1, MA_WRICOMB);

  dims_data1[1] = nvalid;
  this->d_intensities = new CarmaObj<float>(context, dims_data1);

  this->d_fluxPerSub = new CarmaObj<float>(context, dims_data1);
  dims_data1[1] = npix;
  this->d_validsubsx = new CarmaObj<int>(context, dims_data1);
  this->d_validsubsy = new CarmaObj<int>(context, dims_data1);
  dims_data2[1] = nfft;
  dims_data2[2] = nfft;
  this->d_modu_gather = new CarmaObj<float>(current_context, dims_data2);

}

SutraWfs_PyrHR::SutraWfs_PyrHR(
    CarmaContext *context, SutraTelescope *d_tel,
    CarmaObj<cuFloatComplex> *d_camplipup,
    CarmaObj<cuFloatComplex> *d_camplifoc,
    CarmaObj<cuFloatComplex> *d_fttotim, long nxsub, long nvalid, long npupils,
    long npix, long nphase, long nrebin, long nfft, long ntot, long npup,
    float pdiam, float nphotons, float nphot4imat, int lgs, bool fakecam,
    int max_flux_per_pix, int max_pix_value, bool roket, int nbdevices, int *devices)
    : SutraWfs_PyrHR(context, d_tel, d_camplipup, d_camplifoc, d_fttotim,
                          nxsub, nvalid, npupils, npix, nphase, nrebin, nfft,
                          ntot, npup, pdiam, nphotons, nphot4imat, lgs, fakecam,
                          max_flux_per_pix, max_pix_value, roket, devices[0]) {
  long dims_data2[3];
  dims_data2[0] = 2;

  d_hrimg_ngpu.push_back(this->d_hrimg);
  d_camplipup_ngpu.push_back(this->d_camplipup);
  d_camplifoc_ngpu.push_back(this->d_camplifoc);
  d_pyrfocalplane_ngpu.push_back(this->d_pyrfocalplane);
  d_phalfxy_ngpu.push_back(this->d_phalfxy);
  d_fttotim_ngpu.push_back(this->d_fttotim);
  d_screen_ngpu.push_back(
      nullptr);  // init in the SutraWfs_PyrHR::wfs_initarrays
  d_pupil_ngpu.push_back(this->d_pupil);
  d_submask_ngpu.push_back(this->d_submask);
  dims_data2[1] = nfft;
  dims_data2[2] = nfft;
  this->d_modu_gather = new CarmaObj<float>(current_context, dims_data2);
  for (int device = 1; device < nbdevices; device++) {
    current_context->set_active_device(device, 1);
    dims_data2[1] = nfft;
    dims_data2[2] = nfft;
    d_hrimg_ngpu.push_back(
        new CarmaObj<float>(context, dims_data2));  // Useless for SH
    d_camplipup_ngpu.push_back(
        new CarmaObj<cuFloatComplex>(context, dims_data2));
    d_camplifoc_ngpu.push_back(
        new CarmaObj<cuFloatComplex>(context, dims_data2));
    d_pyrfocalplane_ngpu.push_back(new CarmaObj<float>(context, dims_data2));
    cufftHandle *plan =
        this->d_camplipup_ngpu[device]->get_plan();  ///< FFT plan
    carmafft_safe_call(
        cufftPlan2d(plan, dims_data2[1], dims_data2[2], CUFFT_C2C));
    d_phalfxy_ngpu.push_back(
        new CarmaObj<cuFloatComplex>(context, dims_data2));
    d_submask_ngpu.push_back(new CarmaObj<float>(context, dims_data2));
    d_fttotim_ngpu.push_back(
        new CarmaObj<cuFloatComplex>(context, dims_data2));

    dims_data2[1] = ntot;
    dims_data2[2] = ntot;
    d_pupil_ngpu.push_back(new CarmaObj<float>(context, dims_data2));
    d_screen_ngpu.push_back(new CarmaObj<float>(context, dims_data2));

  }
}

SutraWfs_PyrHR::~SutraWfs_PyrHR() {
  for (std::vector<CarmaObj<cuFloatComplex> *>::iterator it =
           this->d_camplipup_ngpu.begin();
       this->d_camplipup_ngpu.end() != it; ++it) {
    if (*it != this->d_camplipup) {
      current_context->set_active_device((*it)->get_device(), 1);
      delete *it;
    }
  }
  this->d_camplipup_ngpu.clear();

  for (std::vector<CarmaObj<cuFloatComplex> *>::iterator it =
           this->d_camplifoc_ngpu.begin();
       this->d_camplifoc_ngpu.end() != it; ++it) {
    if (*it != this->d_camplifoc) {
      current_context->set_active_device((*it)->get_device(), 1);
      delete *it;
    }
  }
  this->d_camplifoc_ngpu.clear();

  for (std::vector<CarmaObj<float> *>::iterator it =
           this->d_pyrfocalplane_ngpu.begin();
       this->d_pyrfocalplane_ngpu.end() != it; ++it) {
    if (*it != this->d_pyrfocalplane) {
      current_context->set_active_device((*it)->get_device(), 1);
      delete *it;
    }
  }
  this->d_pyrfocalplane_ngpu.clear();

  for (std::vector<CarmaObj<cuFloatComplex> *>::iterator it =
           this->d_phalfxy_ngpu.begin();
       this->d_phalfxy_ngpu.end() != it; ++it) {
    if (*it != this->d_phalfxy) {
      current_context->set_active_device((*it)->get_device(), 1);
      delete *it;
    }
  }
  this->d_phalfxy_ngpu.clear();

  for (std::vector<CarmaObj<cuFloatComplex> *>::iterator it =
           this->d_fttotim_ngpu.begin();
       this->d_fttotim_ngpu.end() != it; ++it) {
    if (*it != this->d_fttotim) {
      current_context->set_active_device((*it)->get_device(), 1);
      delete *it;
    }
  }
  this->d_fttotim_ngpu.clear();

  for (std::vector<CarmaObj<float> *>::iterator it =
           this->d_pupil_ngpu.begin();
       this->d_pupil_ngpu.end() != it; ++it) {
    if (*it != this->d_pupil) {
      current_context->set_active_device((*it)->get_device(), 1);
      delete *it;
    }
  }
  this->d_pupil_ngpu.clear();

  for (std::vector<CarmaObj<float> *>::iterator it =
           this->d_submask_ngpu.begin();
       this->d_submask_ngpu.end() != it; ++it) {
    if (*it != this->d_submask) {
      current_context->set_active_device((*it)->get_device(), 1);
      delete *it;
    }
  }
  this->d_submask_ngpu.clear();

  for (std::vector<CarmaObj<float> *>::iterator it =
           this->d_screen_ngpu.begin();
       this->d_screen_ngpu.end() != it; ++it) {
    if (*it != this->d_gs->d_phase->d_screen) {
      current_context->set_active_device((*it)->get_device(), 1);
      delete *it;
    }
  }
  this->d_screen_ngpu.clear();

  for (std::vector<CarmaObj<float> *>::iterator it =
           this->d_hrimg_ngpu.begin();
       this->d_hrimg_ngpu.end() != it; ++it) {
    if (*it != this->d_hrimg) {
      current_context->set_active_device((*it)->get_device(), 1);
      delete *it;
    }
  }
  this->d_hrimg_ngpu.clear();

  current_context->set_active_device(device, 1);
  if (this->d_camplipup != 0L) delete this->d_camplipup;
  if (this->d_camplifoc != 0L) delete this->d_camplifoc;
  if (this->d_pyrfocalplane != 0L) delete this->d_pyrfocalplane;

  if (this->d_fttotim != 0L) delete this->d_fttotim;

  if (this->d_ftkernel != 0L) delete this->d_ftkernel;

  if (this->d_hrimg != 0L) delete this->d_hrimg;
  if (this->d_bincube != 0L) delete this->d_bincube;
  if (this->d_binimg != 0L) delete this->d_binimg;
  if (this->d_intensities != 0L) delete this->d_intensities;
  if (this->d_offsets != 0L) delete this->d_offsets;
  if (this->d_fluxPerSub != 0L) delete this->d_fluxPerSub;
  if (this->d_sincar != 0L) delete this->d_sincar;
  if (this->d_submask != 0L) delete this->d_submask;
  if (this->d_hrmap != 0L) delete this->d_hrmap;

  if (this->d_slopes != 0L) delete this->d_slopes;

  if (this->image_telemetry != 0L) delete this->image_telemetry;

  if (this->d_phasemap != 0L) delete this->d_phasemap;
  if (this->d_ttprojmat != 0L) delete this->d_ttprojmat;
  if (this->d_ttprojvec != 0L) delete this->d_ttprojvec;
  if (this->d_validsubsx != 0L) delete this->d_validsubsx;
  if (this->d_validsubsy != 0L) delete this->d_validsubsy;

  if (this->d_psum != 0L) delete this->d_psum;
  if (this->d_phalfxy != 0L) delete this->d_phalfxy;
  if (this->d_poffsets != 0L) delete this->d_poffsets;
  if (this->pyr_cx != 0L) delete this->pyr_cx;
  if (this->pyr_cy != 0L) delete this->pyr_cy;
  if (this->pyr_mod_weights != 0L) delete this->pyr_mod_weights;

  if (this->lgs) delete this->d_gs->d_lgs;

  if (this->d_gs != 0L) delete this->d_gs;

  delete this->streams;
}

int SutraWfs_PyrHR::load_arrays(cuFloatComplex *halfxy, float *cx,
                                    float *cy, float *weights, float *sincar,
                                    float *submask, int *validsubsx,
                                    int *validsubsy, int *phasemap,
                                    float *fluxPerSub, float *ttprojmat) {
  for (std::vector<CarmaObj<cuFloatComplex> *>::iterator it =
           this->d_phalfxy_ngpu.begin();
       this->d_phalfxy_ngpu.end() != it; ++it) {
    if (*it != this->d_phalfxy) {
      current_context->set_active_device((*it)->get_device(), 1);
      (*it)->host2device(halfxy);
    }
  }
  for (std::vector<CarmaObj<float> *>::iterator it =
           this->d_submask_ngpu.begin();
       this->d_submask_ngpu.end() != it; ++it) {
    if (*it != this->d_submask) {
      current_context->set_active_device((*it)->get_device(), 1);
      (*it)->host2device(submask);
    }
  }
  for (std::vector<CarmaObj<float> *>::iterator it =
           this->d_pupil_ngpu.begin();
       this->d_pupil_ngpu.end() != it; ++it) {
    if (*it != this->d_pupil) {
      current_context->set_active_device((*it)->get_device(), 1);
      (*it)->copy_from(d_pupil->get_data(), d_pupil->get_nb_elements());
    }
  }
  if (d_screen_ngpu.size() > 0)
    d_screen_ngpu[0] = this->d_gs->d_phase->d_screen;

  current_context->set_active_device(device, 1);
  this->d_phalfxy->host2device(halfxy);
  this->d_submask->host2device(submask);
  this->pyr_cx->fill_from(cx);
  this->pyr_cy->fill_from(cy);
  this->pyr_mod_weights->fill_from(weights);
  this->d_sincar->host2device(sincar);
  this->d_validsubsx->host2device(validsubsx);
  this->d_validsubsy->host2device(validsubsy);
  this->d_phasemap->host2device(phasemap);
  this->d_ttprojmat->host2device(ttprojmat);
  this->d_fluxPerSub->host2device(fluxPerSub);

  return EXIT_SUCCESS;
}

int SutraWfs_PyrHR::set_submask(float *submask) {
  int ngpu = d_screen_ngpu.size();
  if (ngpu < 2) {
    this->d_submask->host2device(submask);
  } else {
    for (std::vector<CarmaObj<float> *>::iterator it =
             this->d_submask_ngpu.begin();
         this->d_submask_ngpu.end() != it; ++it) {
      if (*it != this->d_submask) {
        current_context->set_active_device((*it)->get_device(), 1);
        (*it)->host2device(submask);
      }
    }
    this->d_submask->host2device(submask);
  }
  return EXIT_SUCCESS;
}

void SutraWfs_PyrHR::comp_modulation(int cpt) {
  // int cpt = 0;
  int ngpu = d_screen_ngpu.size();
  if (ngpu < 2) {
    carma_safe_call(
        cudaMemset(this->d_camplipup->get_data(), 0,
                   sizeof(cuFloatComplex) * this->d_camplipup->get_nb_elements()));
    pyr_getpup(
        this->d_camplipup->get_data(), this->d_gs->d_phase->d_screen->get_data(),
        this->d_pupil->get_data(), this->ntot, this->nfft, this->d_gs->lambda,
        (this->pyr_cx->get_data())[cpt], (this->pyr_cy->get_data())[cpt],
        this->current_context->get_device(device));
    CarmaFFT(this->d_camplipup->get_data(), this->d_camplifoc->get_data(), -1,
              *this->d_camplipup->get_plan());

    pyr_submask(this->d_camplifoc->get_data(), this->d_submask->get_data(),
                this->nfft, this->current_context->get_device(device));
    // float fact = 1.0f / this->nfft / this->nfft / this->nfft / 2.0;
    float fact = 1.0f * (this->pyr_mod_weights->get_data())[cpt];
    if (compute_pyrfocalplane) {
      abs2(this->d_pyrfocalplane->get_data(), this->d_camplifoc->get_data(),
           this->nfft * this->nfft, fact,
           this->current_context->get_device(device));
    }
    pyr_submaskpyr(this->d_camplifoc->get_data(), this->d_phalfxy->get_data(),
                   this->nfft, this->current_context->get_device(device));
    CarmaFFT(this->d_camplifoc->get_data(), this->d_fttotim->get_data(), 1,
              *this->d_camplipup->get_plan());

    abs2(this->d_hrimg->get_data(), this->d_fttotim->get_data(),
         this->nfft * this->nfft, fact,
         this->current_context->get_device(device));
  } else {
    int cur_device = cpt % ngpu;
    current_context->set_active_device(cur_device, 1);

    carma_safe_call(cudaMemset(
        this->d_camplipup_ngpu[cur_device]->get_data(), 0,
        2 * sizeof(float) * this->d_camplipup_ngpu[cur_device]->get_nb_elements()));
    pyr_getpup(this->d_camplipup_ngpu[cur_device]->get_data(),
               this->d_screen_ngpu[cur_device]->get_data(),
               this->d_pupil_ngpu[cur_device]->get_data(), this->ntot,
               this->nfft, this->d_gs->lambda, (this->pyr_cx->get_data())[cpt],
               (this->pyr_cy->get_data())[cpt],
               this->current_context->get_device(cur_device));
    CarmaFFT(this->d_camplipup_ngpu[cur_device]->get_data(),
              this->d_camplifoc_ngpu[cur_device]->get_data(), -1,
              *this->d_camplipup_ngpu[cur_device]->get_plan());

    pyr_submask(this->d_camplifoc_ngpu[cur_device]->get_data(),
                this->d_submask_ngpu[cur_device]->get_data(), this->nfft,
                this->current_context->get_device(cur_device));
    // float fact = 1.0f / this->nfft / this->nfft / this->nfft / 2.0;
    float fact = 1.0f * (this->pyr_mod_weights->get_data())[cpt];
    if (compute_pyrfocalplane) {
      abs2(this->d_pyrfocalplane_ngpu[cur_device]->get_data(),
           this->d_camplifoc_ngpu[cur_device]->get_data(),
           this->nfft * this->nfft, fact,
           this->current_context->get_device(cur_device));
    }

    pyr_submaskpyr(this->d_camplifoc_ngpu[cur_device]->get_data(),
                   this->d_phalfxy_ngpu[cur_device]->get_data(), this->nfft,
                   this->current_context->get_device(cur_device));
    CarmaFFT(this->d_camplifoc_ngpu[cur_device]->get_data(),
              this->d_fttotim_ngpu[cur_device]->get_data(), 1,
              *this->d_camplipup_ngpu[cur_device]->get_plan());
    abs2(this->d_hrimg_ngpu[cur_device]->get_data(),
         this->d_fttotim_ngpu[cur_device]->get_data(), this->nfft * this->nfft,
         fact, this->current_context->get_device(cur_device));
  }
}

//////////////////////////////
// PYRAMID WAVEFRONT SENSOR //
//////////////////////////////

// It starts by looking for the type of sensor. By default it assumes
// a pyramid wfs. The pyramid can also be explicitely asked for, or
// a roof prism can be asked for as well.
int SutraWfs_PyrHR::comp_generic() {
  /*
   //___________________________________________________________________
   //  PYRAMID SENSOR MODEL

   This code generates pupil images as seen from behind a pyramid wavefront
   sensor algorithm: for (i=0;i<mod_pts;i++) { get phase and multiply by
   exp(i*modu) do fft apply field stop myltiply by exp(i*pyramid) where pyramid
   is the pyramid shape do fft-1 do abs2 add to previous modulation image
   }
   do fft
   multiply by sinc (pixes transfer function)
   take 1 pixels over nrebin pixels in the image
   normalize
   add noise
   */
  current_context->set_active_device(device, 1);
  carma_safe_call(cudaMemset(this->d_binimg->get_data(), 0,
                           sizeof(float) * this->d_binimg->get_nb_elements()));
  for (std::vector<CarmaObj<float> *>::iterator it =
           this->d_screen_ngpu.begin();
       this->d_screen_ngpu.end() != it; ++it) {
    CarmaObj<float> *tmp_screen = this->d_gs->d_phase->d_screen;
    if (*it != tmp_screen) {
      current_context->set_active_device((*it)->get_device(), 1);
      (*it)->copy_from(tmp_screen->get_data(), tmp_screen->get_nb_elements());
    }
  }

  current_context->set_active_device(device, 1);
  carma_safe_call(cudaMemset(this->d_hrimg->get_data(), 0,
                           sizeof(float) * this->d_hrimg->get_nb_elements()));
  for (std::vector<CarmaObj<float> *>::iterator it =
           this->d_hrimg_ngpu.begin();
       this->d_hrimg_ngpu.end() != it; ++it) {
    if (*it != d_hrimg) {
      current_context->set_active_device((*it)->get_device(), 1);
      carma_safe_call(
          cudaMemset((*it)->get_data(), 0, sizeof(float) * (*it)->get_nb_elements()));
    }
  }

  if (compute_pyrfocalplane) {
    current_context->set_active_device(device, 1);
    carma_safe_call(
        cudaMemset(this->d_pyrfocalplane->get_data(), 0,
                   sizeof(float) * this->d_pyrfocalplane->get_nb_elements()));
    for (std::vector<CarmaObj<float> *>::iterator it =
             this->d_pyrfocalplane_ngpu.begin();
         this->d_pyrfocalplane_ngpu.end() != it; ++it) {
      if (*it != d_pyrfocalplane) {
        current_context->set_active_device((*it)->get_device(), 1);
        carma_safe_call(cudaMemset((*it)->get_data(), 0,
                                 sizeof(float) * (*it)->get_nb_elements()));
      }
    }
  }

  current_context->set_active_device(device, 1);

  for (int cpt = 0; cpt < this->npup; cpt++) {
    comp_modulation(cpt);
  }

  current_context->set_active_device(device, 1);

  for (std::vector<CarmaObj<float> *>::iterator it =
           this->d_hrimg_ngpu.begin();
       this->d_hrimg_ngpu.end() != it; ++it) {
    if (*it != d_hrimg) {
      current_context->set_active_device((*it)->get_device(), 1);
      cudaStreamSynchronize(0);
      current_context->set_active_device(device, 1);
      // if (current_context->can_p2p(d_hrimg->get_device(), (*it)->get_device())) {
      //   d_hrimg->axpy(1.0f, (*it), 1, 1);
      // } else {
        d_modu_gather->copy_from((*it)->get_data(), (*it)->get_nb_elements());
        d_hrimg->axpy(1.0f, d_modu_gather, 1, 1);
      // }
    }
  }

  if (compute_pyrfocalplane) {
    for (std::vector<CarmaObj<float> *>::iterator it =
             this->d_pyrfocalplane_ngpu.begin();
         this->d_pyrfocalplane_ngpu.end() != it; ++it) {
      if (*it != d_pyrfocalplane) {
        // if (current_context->can_p2p(d_pyrfocalplane->get_device(),
        //                             (*it)->get_device())) {
        //   d_pyrfocalplane->axpy(1.0f, (*it), 1, 1);
        // } else {
          d_modu_gather->copy_from((*it)->get_data(), (*it)->get_nb_elements());
          d_pyrfocalplane->axpy(1.0f, d_modu_gather, 1, 1);
        // }
      }
    }
  }

  cfillrealp(this->d_fttotim->get_data(), this->d_hrimg->get_data(),
             this->d_hrimg->get_nb_elements(),
             this->current_context->get_device(device));

  CarmaFFT(this->d_fttotim->get_data(), this->d_fttotim->get_data(), -1,
            *this->d_camplipup->get_plan());

  pyr_submask(this->d_fttotim->get_data(), this->d_sincar->get_data(), this->nfft,
              this->current_context->get_device(device));

  CarmaFFT(this->d_fttotim->get_data(), this->d_fttotim->get_data(), 1,
            *this->d_camplipup->get_plan());

  cgetrealp(this->d_hrimg->get_data(), this->d_fttotim->get_data(),
            this->d_hrimg->get_nb_elements(),
            this->current_context->get_device(device));

  pyr_fillbinimg(this->d_binimg->get_data(), this->d_hrimg->get_data(),
                 this->nfft / this->nrebin, this->nfft, this->nrebin, false,
                 this->current_context->get_device(device));
  /*
        pyr_intensities(this->d_intensities->get_data(),
     this->d_binimg->get_data(), this->d_validsubsx->get_data(),
     this->d_validsubsy->get_data(), this->nfft / this->nrebin, this->nvalid,
                        this->current_context->get_device(device));
  */
  int blocks, threads;
  //  get_num_blocks_and_threads(current_context->get_device(device),
  //  this->d_binimg->get_nb_elements(),
  //      blocks, threads);
  this->current_context->set_active_device(device, 1);
  sum_get_num_blocks_and_threads(this->d_binimg->get_nb_elements(),
                            this->current_context->get_device(device), blocks,
                            threads);

  this->d_psum->reset();
  // DEBUG_TRACE("threads %d blocks %d",threads,blocks);
  // reduce(this->d_binimg->get_nb_elements(), threads, blocks,
  // this->d_binimg->get_data(),
  //        this->d_psum->get_data());

  float p_sum =
      reduce<float>(this->d_binimg->get_data(), this->d_binimg->get_nb_elements());

  pyr_fact(this->d_binimg->get_data(), this->nphot / p_sum,
           this->nfft / this->nrebin, 1,
           this->current_context->get_device(device));

  if (this->roket) {  // Get here the binimg before adding noise, usefull
                      // for error budget
    this->d_binimg->copy_into(this->d_binimg_notnoisy->get_data(),
                             this->d_binimg->get_nb_elements());
  }
  // add noise
  if (this->noise > -1) {
    // cout << "adding poisson noise" << endl;
    this->d_binimg->prng('P');
  }
  if (this->noise > 0) {
    // cout << "adding detector noise" << endl;
    this->d_binimg->prng('N', this->noise, 1.0f);
  }

  // Done in getpyr
  //  pyr_intensities(this->d_intensities->get_data(), this->d_binimg->get_data(),
  //             this->d_validsubsx->get_data(), this->d_validsubsy->get_data(),
  //             this->nfft / this->nrebin, this->nvalid,
  //             this->current_context->get_device(device));

  return EXIT_SUCCESS;
}

int SutraWfs_PyrHR::comp_image(bool noise) {
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
  if (this->fakecam)
    result *= digitalize(this->d_camimg->get_data(), this->d_binimg->get_data(),
                         this->d_dark->get_data(), this->d_flat->get_data(),
                         this->max_flux_per_pix, this->max_pix_value,
                         this->d_binimg->get_nb_elements(),
                         this->current_context->get_device(this->device));
  // result *= this->fill_binimage();
  return result;
}

int SutraWfs_PyrHR::fill_binimage(int async) {
  if (this->d_binimg == NULL) {
    DEBUG_TRACE(
        "ERROR : d_bincube not initialized, did you do the allocate_buffers?");
    throw "ERROR : d_bincube not initialized, did you do the allocate_buffers?";
  }
  if (noise > 0) this->d_binimg->prng('N', this->noise);

  this->current_context->set_active_device(device, 1);
  if (async) {
    //    fillbinimg_async(this->image_telemetry, this->d_binimg->get_data(),
    //        this->d_bincube->get_data(), this->npix, this->nvalid_tot,
    //        this->npix * this->nxsub, this->d_validsubsx->get_data(),
    //        this->d_validsubsy->get_data(), this->d_binimg->get_nb_elements(), false,
    //        this->current_context->get_device(device));
    DEBUG_TRACE("ERROR : async version of fill_binimage not implemented...");
    throw "ERROR : async version of fill_binimage not implemented...";
  } else {
    pyr_fillbinimg(this->d_binimg->get_data(), this->d_bincube->get_data(),
                   this->nfft / this->nrebin, false,
                   this->current_context->get_device(device));
  }

  return EXIT_SUCCESS;
}

int SutraWfs_PyrHR::copy_valid_pix(float *img, int *validx, int *validy,
                                      int im_dim) {
  current_context->set_active_device(device, 1);
  copy_imgin_binimg(this->d_binimg->get_data(), this->d_validsubsx->get_data(),
                  this->d_validsubsy->get_data(), this->d_binimg->get_dims(1),
                  img, validx, validy, im_dim, this->d_validsubsx->get_dims(1),
                  this->current_context->get_device(device));
  return EXIT_SUCCESS;
}

int SutraWfs_PyrHR::set_pyr_mod_weights(float *weights, int npts) {
  if (this->npup != npts) {
    DEBUG_TRACE("Number of elements mismatch the modulation points one");
    return EXIT_FAILURE;
  }
  this->pyr_mod_weights->fill_from(weights);
  return EXIT_SUCCESS;
}

int SutraWfs_PyrHR::set_pyr_modulation_points(float *cx, float *cy,
                                            float *weights, int npts) {
  int status;
  status = this->set_pyr_modulation_points(cx, cy, npts);
  status *= this->set_pyr_mod_weights(weights, npts);
  return status;
}

int SutraWfs_PyrHR::set_pyr_modulation_points(float *cx, float *cy, int npts) {
  current_context->set_active_device(device, 1);
  this->npup = npts;
  if (this->pyr_cx != 0L) {
    delete this->pyr_cx;
    delete this->pyr_cy;
  }
  if (this->pyr_mod_weights != 0L) {
    delete this->pyr_mod_weights;
  }
  long dims_data1[2] = {1, npts};
  this->pyr_cx = new CarmaHostObj<float>(dims_data1, cx, MA_WRICOMB);
  this->pyr_cy = new CarmaHostObj<float>(dims_data1, cy, MA_WRICOMB);
  this->pyr_mod_weights = new CarmaHostObj<float>(dims_data1, MA_WRICOMB);
  this->pyr_mod_weights->fill(1.f);

  return EXIT_SUCCESS;
}

int SutraWfs_PyrHR::comp_nphot(float ittime, float optthroughput,
                                    float diam, float cobs, float zerop,
                                    float gsmag) {
  this->d_gs->mag = gsmag;
  this->nphot = zerop * pow(10., -0.4 * gsmag) * ittime * optthroughput *
                CARMA_PI / 4. * (1 - cobs * cobs) * diam * diam;
  return EXIT_SUCCESS;
}

int SutraWfs_PyrHR::set_phalfxy(cuFloatComplex *phalfxy) {
  int ngpus = this->d_screen_ngpu.size();
  for (int dev = 0; dev < ngpus; dev++) {
      current_context->set_active_device(dev, 1);
      this->d_phalfxy_ngpu[dev]->host2device(phalfxy);
  }
  current_context->set_active_device(this->device, 1);
  return EXIT_SUCCESS;
}
