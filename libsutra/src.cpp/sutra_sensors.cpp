// This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
//
// COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU
// Lesser General Public License as published by the Free Software Foundation, either version 3 of
// the License, or any later version.
//
// COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
// even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License along with COMPASS. If
// not, see <https://www.gnu.org/licenses/>
//
//  Copyright (C) 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>

//! \file      sutra_sensors.cpp
//! \ingroup   libsutra
//! \class     SutraSensors
//! \brief     this class provides the sensors features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#include <sutra_sensors.hpp>

SutraSensors::SutraSensors(CarmaContext *context, SutraTelescope *d_tel,
                             vector<string> type, int32_t nwfs, int64_t *nxsub,
                             int64_t *nvalid, int64_t *npupils, int64_t *npix,
                             int64_t *nphase, int64_t *nrebin, int64_t *nfft, int64_t *ntot,
                             int64_t *npup, float *pdiam, float *nphot,
                             float *nphot4imat, int32_t *lgs, bool *fakecam,
                             int32_t *max_flux_per_pix, int32_t *max_pix_value, int32_t device,
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
    int32_t maxnfft = nfft[0];
    int32_t maxntot = ntot[0];
    int32_t maxnvalid = 0;
    int32_t maxnvalid_tot = nvalid[0];
    int32_t is_lgs = (lgs[0] > 0 ? 1 : 0);
    for (int32_t i = 1; i < nwfs; i++) {
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
    int64_t dims_data3[4] = {
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
  for (int32_t i = 0; i < nwfs; i++) {
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
      const int32_t ngpu = context->get_ndevice();
      // DEBUG_TRACE("using pyrhr with %d GPUs", ngpu);
      if (ngpu == 1) {
        wfs = new SutraWfs_PyrHR(
            context, d_tel, this->d_camplipup, this->d_camplifoc,
            this->d_fttotim, nxsub[i], nvalid[i], npupils[i], npix[i],
            nphase[i], nrebin[i], nfft[i], ntot[i], npup[i], pdiam[i], nphot[i],
            nphot4imat[i], lgs[i], fakecam[i], max_flux_per_pix[i], max_pix_value[i],
            roket, device);
      } else {
        int32_t devices[ngpu];
        for (int32_t i = 0; i < ngpu; i++) {
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
  map<vector<int32_t>, cufftHandle *>::iterator it;
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

int32_t SutraSensors::allocate_buffers() {
  current_context->set_active_device(device, 1);
  for (size_t idx = 0; idx < this->d_wfs.size(); idx++) {
    this->d_wfs[idx]->allocate_buffers(this->campli_plans, this->fttotim_plans);
  }
  return EXIT_SUCCESS;
}

int32_t SutraSensors::define_mpi_rank(int32_t rank, int32_t size) {
  current_context->set_active_device(device, 1);
  for (size_t idx = 0; idx < this->d_wfs.size(); idx++) {
    this->d_wfs[idx]->define_mpi_rank(rank, size);
  }
  return EXIT_SUCCESS;
}

int32_t SutraSensors::initgs(float *xpos, float *ypos, float *lambda, float *mag,
                          float zerop, int64_t *size, float *noise, int64_t *seed,
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
int32_t SutraSensors::initgs(float *xpos, float *ypos, float *lambda, float *mag,
                          float zerop, int64_t *size, float *noise, float *G,
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
int32_t SutraSensors::initgs(float *xpos, float *ypos, float *lambda, float *mag,
                          float zerop, int64_t *size, float *G, float *thetaML,
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

int32_t SutraSensors::set_noise(int32_t nwfs, float noise, int64_t seed) {
  this->d_wfs[nwfs]->set_noise(noise, seed);
  return EXIT_SUCCESS;
}

int32_t SutraSensors::set_field_stop(int32_t nwfs, float* field_stop, int32_t N) {
  if(this->d_wfs[nwfs]->type == "sh") {
    SutraWfsSH *wfs = dynamic_cast<SutraWfsSH *>(this->d_wfs[nwfs]);
    wfs->set_field_stop(this->campli_plans, field_stop, N);
  }
  return EXIT_SUCCESS;
}