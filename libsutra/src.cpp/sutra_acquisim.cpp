/*
 * sutra_acquisim.cpp
 *
 *  Created on: Oct 28, 2013
 *      Author: sevin
 */

#include "sutra_acquisim.h"

sutra_acquisim::sutra_acquisim(sutra_sensors *sensors, int wfs_num)
    : device(sensors->d_wfs[wfs_num]->device),
      type(sensors->d_wfs[wfs_num]->type),
      nxsub(sensors->d_wfs[wfs_num]->nxsub),
      nvalid(sensors->d_wfs[wfs_num]->nvalid),
      npix(sensors->d_wfs[wfs_num]->npix),
      d_validsubsx(sensors->d_wfs[wfs_num]->d_validsubsx),
      d_validsubsy(sensors->d_wfs[wfs_num]->d_validsubsy) {
  // TODO Auto-generated constructor stub
  if (sensors->d_wfs[wfs_num]->type == "sh")
    this->wfs = dynamic_cast<sutra_wfs_sh *>(sensors->d_wfs[wfs_num]);
}

sutra_acquisim::~sutra_acquisim() {
  // TODO Auto-generated destructor stub
}

int sutra_acquisim::set_validsubs(int64_t nvalid, int32_t *validsubsx,
                                  int32_t *validsubsy) {
  int64_t dims[2] = {1, nvalid};
  this->d_validsubsx =
      new carma_obj<int32_t>(this->current_context, dims, validsubsx);
  this->d_validsubsy =
      new carma_obj<int32_t>(this->current_context, dims, validsubsy);
  return EXIT_SUCCESS;
}

int sutra_acquisim::comp_image(long *dims, float *bimage) {
  return comp_image(dims, bimage, this->wfs->d_bincube);
}

int sutra_acquisim::comp_image(long *dims, float *bimage,
                               carma_obj<float> *d_bincube) {
  this->current_context->set_activeDevice(this->device, 1);
  carma_obj<float> tmp_yObj(this->current_context, dims, bimage);
  return fillbincube<float>(tmp_yObj, *d_bincube, this->npix, this->nvalid,
                            this->npix * this->nxsub, *(this->d_validsubsx),
                            *(this->d_validsubsy),
                            this->current_context->get_device(this->device));
}

int sutra_acquisim::comp_image_2D(long *dims, float *bimage, int *num_ssp) {
  this->current_context->set_activeDevice(this->device, 1);
  carma_obj<float> tmp_yObj(this->current_context, dims, bimage);
  long dims1[2] = {1, this->nxsub * this->nxsub};
  carma_obj<int> d_num_ssp(this->current_context, dims1, num_ssp);
  return fillbincube_2D<float>(tmp_yObj, *(this->wfs->d_bincube), this->npix,
                               this->nxsub, d_num_ssp);
}

int sutra_acquisim::comp_image_tele(long *dims, float *bimage) {
  this->current_context->set_activeDevice(this->device, 1);
  carma_obj<float> tmp_yObj(this->current_context, dims, bimage);
  return fillbincube_async<float>(
      this->wfs->image_telemetry, tmp_yObj, *(this->wfs->d_bincube), this->npix,
      this->nvalid, this->npix * this->nxsub, *(this->d_validsubsx),
      *(this->d_validsubsy), this->wfs->d_binimg->getNbElem(),
      this->current_context->get_device(this->device));
}
