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

//! \file      sutra_sensors.cpp
//! \ingroup   libsutra
//! \class     SutraSensors
//! \brief     this class provides the sensors features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.0.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#include <sutra_sensors.h>

SutraSensors::SutraSensors(CarmaContext *context, SutraTelescope *d_tel,
                             vector<string> type, int nwfs, long *nxsub,
                             long *nvalid, long *npupils, long *npix,
                             long *nphase, long *nrebin, long *nfft, long *ntot,
                             long *npup, float *pdiam, float *nphot,
                             float *nphot4imat, int *lgs, bool *fakecam,
                             int *max_flux_per_pix, int *max_pix_value, int device,
                             bool roket)
    : device(device),
      current_context(context),
      d_camplipup(nullptr),
      d_camplifoc(nullptr),
      d_fttotim(nullptr),
      d_ftlgskern(nullptr),
      d_lgskern(nullptr) {
  current_context->set_active_device(device, 1);
  // DEBUG_TRACE("Before create sensors : ");print_mem_info();
  if (type[0].compare("sh") == 0) {
    int maxnfft = nfft[0];
    int maxntot = ntot[0];
    int maxnvalid = 0;
    int maxnvalid_tot = nvalid[0];
    int is_lgs = (lgs[0] > 0 ? 1 : 0);
    for (int i = 1; i < nwfs; i++) {
      if (type[i].compare("shlo") != 0) {
        if (nvalid[i] > maxnvalid_tot) maxnvalid_tot = nvalid[i];
        if (ntot[i] > maxntot) {
          maxntot = ntot[i];
        }
        if (nfft[i] > maxnfft) {
          maxnfft = nfft[i];
        }
        if (ntot[i] == nfft[i]) {
          if (nvalid[i] > maxnvalid) {
            maxnvalid = nvalid[i];
          }
        }
      }
      if (lgs[i] > 0) is_lgs = 1;
    }
    // DEBUG_TRACE("maxntot : %d maxnfft : %d maxnvalid : %d nmaxhr : %d \n
    // ",maxntot,maxnfft,maxnvalid,compute_nmaxhr(nvalid[wfs4ntot]));
    long dims_data3[4] = {
        3, maxnfft, maxnfft, maxnvalid_tot /*nvalid[wfs4nfft]*/
    };
    this->d_camplifoc = new CarmaObj<cuFloatComplex>(context, dims_data3);
    this->d_camplipup = new CarmaObj<cuFloatComplex>(context, dims_data3);

    dims_data3[1] = maxntot;
    dims_data3[2] = maxntot;
    dims_data3[3] = compute_nmaxhr(maxnvalid_tot /*nvalid[wfs4ntot]*/);
    if (maxnvalid > dims_data3[3]) dims_data3[3] = maxnvalid;
    this->d_fttotim = new CarmaObj<cuFloatComplex>(context, dims_data3);
    if (is_lgs) {
      this->d_ftlgskern = new CarmaObj<cuFloatComplex>(context, dims_data3);
      this->d_lgskern = new CarmaObj<float>(context, dims_data3);
    } else {
      this->d_ftlgskern = nullptr;
      this->d_lgskern = nullptr;
    }
  } else {  // Pyramid case : shared arrays are not used
    this->d_camplifoc = nullptr;
    this->d_camplipup = nullptr;
    this->d_fttotim = nullptr;
    this->d_ftlgskern = nullptr;
    this->d_lgskern = nullptr;
  }

  // DEBUG_TRACE("After creating sensors arrays : ");print_mem_info();
  for (int i = 0; i < nwfs; i++) {
    SutraWfs *wfs = NULL;
    if (type[i].compare("sh") == 0 || type[i].compare("shlo") == 0) {
      bool is_low_order = false;
      if (type[i].compare("shlo") == 0) is_low_order = true;
      wfs = new SutraWfsSH(
          context, d_tel, this->d_camplipup, this->d_camplifoc, this->d_fttotim,
          nxsub[i], nvalid[i], npix[i], nphase[i], nrebin[i], nfft[i], ntot[i],
          npup[i], pdiam[i], nphot[i], nphot4imat[i], lgs[i], fakecam[i],
          max_flux_per_pix[i], max_pix_value[i], is_low_order, roket, device);
    }
    if (type[i].compare("pyrhr") == 0) {
      const int ngpu = context->get_ndevice();
      DEBUG_TRACE("using pyrhr with %d GPUs", ngpu);
      if (ngpu == 1) {
        wfs = new SutraWfs_PyrHR(
            context, d_tel, this->d_camplipup, this->d_camplifoc,
            this->d_fttotim, nxsub[i], nvalid[i], npupils[i], npix[i],
            nphase[i], nrebin[i], nfft[i], ntot[i], npup[i], pdiam[i], nphot[i],
            nphot4imat[i], lgs[i], fakecam[i], max_flux_per_pix[i], max_pix_value[i],
            roket, device);
      } else {
        int devices[ngpu];
        for (int i = 0; i < ngpu; i++) {
          devices[i] = i;
        }
        wfs = new SutraWfs_PyrHR(
            context, d_tel, this->d_camplipup, this->d_camplifoc,
            this->d_fttotim, nxsub[i], nvalid[i], npupils[i], npix[i],
            nphase[i], nrebin[i], nfft[i], ntot[i], npup[i], pdiam[i], nphot[i],
            nphot4imat[i], lgs[i], fakecam[i], max_flux_per_pix[i], max_pix_value[i],
            roket, ngpu, devices);
      }
    }

    d_wfs.push_back(wfs);

    //		DEBUG_TRACE("After creating wfs #%d : ",i);print_mem_info();
  }
  //	DEBUG_TRACE("Final sensors : ");print_mem_info();
  this->allocate_buffers();
}

SutraSensors::~SutraSensors() {
  current_context->set_active_device(device, 1);
  //  for (size_t idx = 0; idx < (this->d_wfs).size(); idx++) {
  while ((this->d_wfs).size() > 0) {
    delete this->d_wfs.back();
    d_wfs.pop_back();
  }
  map<vector<int>, cufftHandle *>::iterator it;
  for (it = campli_plans.begin(); it != campli_plans.end(); it++) {
    cufftDestroy(*it->second);
    free(it->second);
  }
  for (it = fttotim_plans.begin(); it != fttotim_plans.end(); it++) {
    cufftDestroy(*it->second);
    free(it->second);
  }
  for (it = ftlgskern_plans.begin(); it != ftlgskern_plans.end(); it++) {
    cufftDestroy(*it->second);
    free(it->second);
  }

  if (this->d_camplifoc != nullptr) delete this->d_camplifoc;
  if (this->d_camplipup != nullptr) delete this->d_camplipup;
  if (this->d_fttotim != nullptr) delete this->d_fttotim;
  if (this->d_lgskern != nullptr) delete this->d_lgskern;
  if (this->d_ftlgskern != nullptr) delete this->d_ftlgskern;
}

int SutraSensors::allocate_buffers() {
  current_context->set_active_device(device, 1);
  for (size_t idx = 0; idx < this->d_wfs.size(); idx++) {
    this->d_wfs[idx]->allocate_buffers(this->campli_plans, this->fttotim_plans);
  }
  return EXIT_SUCCESS;
}

int SutraSensors::define_mpi_rank(int rank, int size) {
  current_context->set_active_device(device, 1);
  for (size_t idx = 0; idx < this->d_wfs.size(); idx++) {
    this->d_wfs[idx]->define_mpi_rank(rank, size);
  }
  return EXIT_SUCCESS;
}

int SutraSensors::initgs(float *xpos, float *ypos, float *lambda, float *mag,
                          float zerop, long *size, float *noise, long *seed,
                          float *G, float *thetaML, float *dx, float *dy) {
  current_context->set_active_device(device, 1);
  for (size_t idx = 0; idx < this->d_wfs.size(); idx++) {
    this->d_wfs[idx]->wfs_initgs(
        this->d_lgskern, this->d_ftlgskern, this->ftlgskern_plans, xpos[idx],
        ypos[idx], lambda[idx], mag[idx], zerop, size[idx], noise[idx],
        seed[idx], G[idx], thetaML[idx], dx[idx], dy[idx]);
  }
  return EXIT_SUCCESS;
}
int SutraSensors::initgs(float *xpos, float *ypos, float *lambda, float *mag,
                          float zerop, long *size, float *noise, float *G,
                          float *thetaML, float *dx, float *dy) {
  current_context->set_active_device(device, 1);
  for (size_t idx = 0; idx < this->d_wfs.size(); idx++) {
    this->d_wfs[idx]->wfs_initgs(
        this->d_lgskern, this->d_ftlgskern, this->ftlgskern_plans, xpos[idx],
        ypos[idx], lambda[idx], mag[idx], zerop, size[idx], noise[idx],
        1234 + idx, G[idx], thetaML[idx], dx[idx], dy[idx]);
  }
  return EXIT_SUCCESS;
}
int SutraSensors::initgs(float *xpos, float *ypos, float *lambda, float *mag,
                          float zerop, long *size, float *G, float *thetaML,
                          float *dx, float *dy) {
  current_context->set_active_device(device, 1);
  for (size_t idx = 0; idx < this->d_wfs.size(); idx++) {
    this->d_wfs[idx]->wfs_initgs(
        this->d_lgskern, this->d_ftlgskern, this->ftlgskern_plans, xpos[idx],
        ypos[idx], lambda[idx], mag[idx], zerop, size[idx], -1, 1234 + idx,
        G[idx], thetaML[idx], dx[idx], dy[idx]);
  }
  return EXIT_SUCCESS;
}

int SutraSensors::set_noise(int nwfs, float noise, long seed) {
  this->d_wfs[nwfs]->set_noise(noise, seed);
  return EXIT_SUCCESS;
}
