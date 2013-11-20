/*
 * sutra_acquisim.h
 *
 *  Created on: Oct 28, 2013
 *      Author: sevin
 */

#ifndef SUTRA_ACQUISIM_H_
#define SUTRA_ACQUISIM_H_

#include <sutra_wfs.h>

class sutra_acquisim {
public:
	sutra_wfs *wfs;

public:
 sutra_acquisim(sutra_sensors *sensors, int wfs_num);
 sutra_acquisim(const sutra_acquisim& acquisim);
 ~sutra_acquisim();

 int comp_image_tele(long *dims, float *bimage);
 int comp_image(long *dims, float *bimage);

private:
};

// General utilities
int fillbincube(float *bimage, float *bcube, int npix, int nsub, int Nsub, int *ivalid, int *jvalid, int device);
int fillbincube_async(carma_streams *streams, carma_obj<float> *bimage, carma_obj<float> *bcube, int npix, int nsub, int Nsub, int *ivalid, int *jvalid, int device);
int fillbincube_async(carma_host_obj<float> *image_telemetry, float *bimage, float *bcube, int npix, int nsub, int Nsub, int *ivalid, int *jvalid,  int nim, int device);

#endif /* SUTRA_ACQUISIM_H_ */
