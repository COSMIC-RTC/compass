#ifndef _SUTRA_CENTROIDER_TCOG_H_
#define _SUTRA_CENTROIDER_TCOG_H_

#include <sutra_centroider.h>

template <class Tin, class T>
class sutra_centroider_tcog : public sutra_centroider<Tin, T> {
 public:
  float threshold;

 public:
  sutra_centroider_tcog(carma_context *context, sutra_wfs *wfs, long nvalid,
                        float offset, float scale, int device);
  sutra_centroider_tcog(const sutra_centroider_tcog &centroider);
  ~sutra_centroider_tcog();

  string get_type();

  int set_threshold(float threshold);

  int get_cog(float *cube, T *intensities, T *centroids, int nvalid, int npix,
              int ntot);
  int get_cog(T *intensities, T *slopes, bool noise);
  int get_cog();
};

template <class T>
void get_centroids(int size, int threads, int blocks, int n, float *d_idata,
                   T *d_odata, T *ref, int *validx, int *validy, T *intensities,
                   float threshold, float scale, float offset,
                   carma_device *device);

#endif  // _SUTRA_CENTROIDER_TCOG_H_
