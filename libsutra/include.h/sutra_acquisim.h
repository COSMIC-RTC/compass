/**
 * \file sutra_acquisim.h
 *
 * \class sutra_acquisim
 *
 * \ingroup libsutra
 *
 * \brief this class provides the acquisition simulator to COMPASS
 *
 * \authors Damien Gratadour & Arnaud Sevin & Florian Ferreira
 *
 * \version 4.2.0
 *
 * \date 2011/01/28
 *
 */

#ifndef SUTRA_ACQUISIM_H_
#define SUTRA_ACQUISIM_H_

#include <sutra_sensors.h>
#include <sutra_wfs_sh.h>

class sutra_acquisim {
 public:
  int device;
  string type;
  long nxsub;
  long nvalid;
  long npix;

  carma_context *current_context;

  sutra_wfs_sh *wfs;
  carma_obj<int32_t> *d_validsubsx;
  carma_obj<int32_t> *d_validsubsy;

 public:
  sutra_acquisim(sutra_sensors *sensors, int wfs_num);
  sutra_acquisim(const sutra_acquisim &acquisim);
  ~sutra_acquisim();

  int set_validsubs(int64_t nvalid, int32_t *validsubsx, int32_t *validsubsy);

  int comp_image_tele(long *dims, float *bimage);
  int comp_image(long *dims, float *bimage);
  int comp_image(long *dims, float *bimage, carma_obj<float> *d_bincube);
  int comp_image_2D(long *dims, float *bimage, int *num_ssp);

 private:
};

// General utilities
template <class T>
int fillbincube_2D(T *bimage, T *bcube, int npix, int nsub, int *valid);

template <class T>
int fillbincube(T *bimage, T *bcube, int npix, int nsub, int Nsub, int *ivalid,
                int *jvalid, carma_device *device);

template <class T>
int fillbincube_async(carma_streams *streams, carma_obj<T> *bimage,
                      carma_obj<T> *bcube, int npix, int nsub, int Nsub,
                      int *ivalid, int *jvalid, carma_device *device);

template <class T>
int fillbincube_async(carma_host_obj<T> *image_telemetry, T *bimage, T *bcube,
                      int npix, int nsub, int Nsub, int *ivalid, int *jvalid,
                      int nim, carma_device *device);

#endif /* SUTRA_ACQUISIM_H_ */
