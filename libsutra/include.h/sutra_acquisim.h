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
  int comp_image_2D(long *dims, float *bimage, int *num_ssp);

private:
};

// General utilities
template<class T>
int fillbincube_2D(T *bimage, T *bcube, int npix, int nsub, int *valid);

template<class T>
int fillbincube(T *bimage, T *bcube, int npix, int nsub, int Nsub, int *ivalid,
    int *jvalid, carma_device *device);

template<class T>
int fillbincube_async(carma_streams *streams, carma_obj<T> *bimage,
    carma_obj<T> *bcube, int npix, int nsub, int Nsub, int *ivalid, int *jvalid,
    carma_device *device);

template<class T>
int fillbincube_async(carma_host_obj<T> *image_telemetry, T *bimage, T *bcube,
    int npix, int nsub, int Nsub, int *ivalid, int *jvalid, int nim,
    carma_device *device);

#endif /* SUTRA_ACQUISIM_H_ */
