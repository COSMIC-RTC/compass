#ifndef _SUTRA_CENTROIDER_MASKEDPIX_H_
#define _SUTRA_CENTROIDER_MASKEDPIX_H_

#include <sutra_centroider.h>
#include <sutra_wfs_pyr_pyrhr.h>
#include <string>

class sutra_centroider_maskedPix : public sutra_centroider {
 public:
  sutra_centroider_maskedPix(carma_context *context, sutra_wfs *wfs,
                             long nvalid, long npupils, float offset,
                             float scale, int device);

  ~sutra_centroider_maskedPix();

  string get_type();

  int get_maskedPix(float *cube, float *intensities, float *centroids,
                    int *subindx, int *subindy, int nvalid, int ns, int nim);
  int get_cog(float *cube, float *intensities, float *centroids, int nvalid,
              int npix, int ntot);
  int get_cog(float *intensities, float *slopes, bool noise);
  int get_cog();
};

template <class T>
void fill_intensities(T *intensities, T *cube, int *subindx, int *subindy,
                      int ns, int nslopes, carma_device *device);
template <class T>
void getMaskedPix(T *centroids, T *img, int *subindx, int *subindy, T *psum,
                  int ns, int nslopes, carma_device *device);

#endif  // _SUTRA_CENTROIDER_MASKEDPIX_H_
