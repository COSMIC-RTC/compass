/*
 * sutra_acquisim.cpp
 *
 *  Created on: Oct 28, 2013
 *      Author: sevin
 */

#include "sutra_acquisim.h"

sutra_acquisim::sutra_acquisim(sutra_sensors *sensors, int wfs_num) {
	// TODO Auto-generated constructor stub
	this->wfs=sensors->d_wfs[wfs_num];
}

sutra_acquisim::~sutra_acquisim() {
	// TODO Auto-generated destructor stub
}

int sutra_acquisim::comp_image(long *dims, float *bimage)
{
	carma_obj<float> tmp_yObj(this->wfs->current_context, dims, bimage);
	return fillbincube(tmp_yObj.getData(), wfs->d_bincube->getData(), wfs->npix , wfs->nvalid, wfs->npix * wfs->nxsub, wfs->d_validsubsx->getData(), wfs->d_validsubsy->getData(), wfs->device);

}

int sutra_acquisim::comp_image_tele(long *dims, float *bimage)
{
	carma_obj<float> tmp_yObj(this->wfs->current_context, dims, bimage);
	return fillbincube_async(wfs->image_telemetry, tmp_yObj.getData(), wfs->d_bincube->getData(), wfs->npix,
			wfs->nvalid, wfs->npix * wfs->nxsub, wfs->d_validsubsx->getData(),
			wfs->d_validsubsy->getData(), wfs->d_binimg->getNbElem(), wfs->device);
}

