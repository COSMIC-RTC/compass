/*
 * sutra_acquisim.cpp
 *
 *  Created on: Oct 28, 2013
 *      Author: sevin
 */

#include "sutra_acquisim.h"

sutra_acquisim::sutra_acquisim(sutra_sensors *sensors, int wfs_num) {
  // TODO Auto-generated constructor stub
  this->wfs = sensors->d_wfs[wfs_num];
}

sutra_acquisim::~sutra_acquisim() {
  // TODO Auto-generated destructor stub
}

int sutra_acquisim::comp_image(long *dims, float *bimage) {
  wfs->current_context->set_activeDevice(wfs->device,1);
  carma_obj<float> tmp_yObj(this->wfs->current_context, dims, bimage);
  return fillbincube<float>(tmp_yObj, *(wfs->d_bincube), wfs->npix, wfs->nvalid,
      wfs->npix * wfs->nxsub, *(wfs->d_validsubsx), *(wfs->d_validsubsy),
      wfs->current_context->get_device(wfs->device));
}

int sutra_acquisim::comp_image_2D(long *dims, float *bimage, int * num_ssp) {
  wfs->current_context->set_activeDevice(wfs->device,1);
  carma_obj<float> tmp_yObj(this->wfs->current_context, dims, bimage);
  long dims1[2] = {1, wfs->nxsub*wfs->nxsub};
  carma_obj<int> d_num_ssp(this->wfs->current_context, dims1, num_ssp);
  return fillbincube_2D<float>(tmp_yObj, *(wfs->d_bincube), wfs->npix, wfs->nxsub, d_num_ssp);
}

int sutra_acquisim::comp_image_tele(long *dims, float *bimage) {
  wfs->current_context->set_activeDevice(wfs->device,1);
  carma_obj<float> tmp_yObj(this->wfs->current_context, dims, bimage);
  return fillbincube_async<float>(wfs->image_telemetry, tmp_yObj,
      *(wfs->d_bincube), wfs->npix, wfs->nvalid, wfs->npix * wfs->nxsub,
      *(wfs->d_validsubsx), *(wfs->d_validsubsy), wfs->d_binimg->getNbElem(),
      wfs->current_context->get_device(wfs->device));
}

