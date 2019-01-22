#ifndef _SUTRA_CENTROIDER_MASKEDPIX_H_
#define _SUTRA_CENTROIDER_MASKEDPIX_H_

#include <sutra_centroider.h>
#include <sutra_wfs_pyr_pyrhr.h>
#include <string>

template <typename T>
class sutra_centroider_maskedPix : public sutra_centroider<T> {
 public:
  sutra_centroider_maskedPix(carma_context *context, sutra_wfs *wfs,
                             long nvalid, long npupils, T offset, T scale,
                             int device);

  ~sutra_centroider_maskedPix();

  string get_type();

  int get_maskedPix(T *img, T *intensities, T *centroids, int *subindx,
                    int *subindy, int nvalid, int ns);
  int get_cog(T *img, T *intensities, T *centroids, int nvalid, int npix,
              int ntot);
  int get_cog(T *intensities, T *slopes, bool noise);
  int get_cog();
};

template <class T>
void fill_intensities(T *intensities, T *img, int *subindx, int *subindy,
                      int ns, int nslopes, carma_device *device);
template <class T>
void getMaskedPix(T *centroids, T *ref, T *img, int *subindx, int *subindy,
                  T *psum, int ns, int nslopes, carma_device *device);

#endif  // _SUTRA_CENTROIDER_MASKEDPIX_H_
