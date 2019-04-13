#ifndef _SUTRA_CENTROIDER_MASKEDPIX_H_
#define _SUTRA_CENTROIDER_MASKEDPIX_H_

#include <sutra_centroider.h>
#include <sutra_wfs_pyr_pyrhr.h>
#include <string>

template <class Tin, class T>
class sutra_centroider_maskedPix : public sutra_centroider<Tin, T> {
 public:
  sutra_centroider_maskedPix(carma_context *context, sutra_wfs *wfs,
                             long nvalid, long npupils, float offset,
                             float scale, bool filter_TT, int device);

  ~sutra_centroider_maskedPix();

  string get_type();

  int get_maskedPix(float *img, float *intensities, T *centroids, int *subindx,
                    int *subindy, int nvalid, int ns);
  int get_cog(float *img, float *intensities, T *centroids, int nvalid,
              int npix, int ntot);
  int get_cog(float *intensities, T *slopes, bool noise);
  int get_cog();
};

void fill_intensities(float *intensities, float *img, int *subindx,
                      int *subindy, int ns, int nslopes, carma_device *device);
template <class T>
void getMaskedPix(T *centroids, T *ref, float *img, int *subindx, int *subindy,
                  float *psum, int ns, int nslopes, carma_device *device);

#endif  // _SUTRA_CENTROIDER_MASKEDPIX_H_
