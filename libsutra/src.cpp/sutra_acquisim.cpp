// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_acquisim.cpp
//! \ingroup   libsutra
//! \class     SutraAcquisim
//! \brief     this class provides the acquisition simulator to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.3
//! \date      2022/01/24

#include "sutra_acquisim.h"

SutraAcquisim::SutraAcquisim(SutraSensors *sensors, int wfs_num)
    : device(sensors->d_wfs[wfs_num]->device),
      type(sensors->d_wfs[wfs_num]->type),
      nxsub(sensors->d_wfs[wfs_num]->nxsub),
      nvalid(sensors->d_wfs[wfs_num]->nvalid),
      npix(sensors->d_wfs[wfs_num]->npix),
      d_validsubsx(sensors->d_wfs[wfs_num]->d_validsubsx),
      d_validsubsy(sensors->d_wfs[wfs_num]->d_validsubsy) {
  // TODO Auto-generated constructor stub
  if (sensors->d_wfs[wfs_num]->type == "sh")
    this->wfs = dynamic_cast<SutraWfsSH *>(sensors->d_wfs[wfs_num]);
}

SutraAcquisim::~SutraAcquisim() {
  // TODO Auto-generated destructor stub
}

int SutraAcquisim::set_validsubs(int64_t nvalid, int32_t *validsubsx,
                                  int32_t *validsubsy) {
  int64_t dims[2] = {1, nvalid};
  this->d_validsubsx =
      new CarmaObj<int32_t>(this->current_context, dims, validsubsx);
  this->d_validsubsy =
      new CarmaObj<int32_t>(this->current_context, dims, validsubsy);
  return EXIT_SUCCESS;
}

int SutraAcquisim::comp_image(long *dims, float *bimage) {
  return comp_image(dims, bimage, this->wfs->d_bincube);
}

int SutraAcquisim::comp_image(long *dims, float *bimage,
                               CarmaObj<float> *d_bincube) {
  this->current_context->set_active_device(this->device, 1);
  CarmaObj<float> tmp_yObj(this->current_context, dims, bimage);
  return fillbincube<float>(tmp_yObj, *d_bincube, this->npix, this->nvalid,
                            this->npix * this->nxsub, *(this->d_validsubsx),
                            *(this->d_validsubsy),
                            this->current_context->get_device(this->device));
}

int SutraAcquisim::comp_image_2D(long *dims, float *bimage, int *num_ssp) {
  this->current_context->set_active_device(this->device, 1);
  CarmaObj<float> tmp_yObj(this->current_context, dims, bimage);
  long dims1[2] = {1, this->nxsub * this->nxsub};
  CarmaObj<int> d_num_ssp(this->current_context, dims1, num_ssp);
  return fillbincube_2D<float>(tmp_yObj, *(this->wfs->d_bincube), this->npix,
                               this->nxsub, d_num_ssp);
}

int SutraAcquisim::comp_image_tele(long *dims, float *bimage) {
  this->current_context->set_active_device(this->device, 1);
  CarmaObj<float> tmp_yObj(this->current_context, dims, bimage);
  return fillbincube_async<float>(
      this->wfs->image_telemetry, tmp_yObj, *(this->wfs->d_bincube), this->npix,
      this->nvalid, this->npix * this->nxsub, *(this->d_validsubsx),
      *(this->d_validsubsy), this->wfs->d_binimg->get_nb_elements(),
      this->current_context->get_device(this->device));
}
