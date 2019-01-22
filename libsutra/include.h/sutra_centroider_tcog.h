#ifndef _SUTRA_CENTROIDER_TCOG_H_
#define _SUTRA_CENTROIDER_TCOG_H_

#include <sutra_centroider.h>

template <typename T>
class sutra_centroider_tcog : public sutra_centroider<T> {
 public:
  T threshold;

 public:
  sutra_centroider_tcog(carma_context *context, sutra_wfs *wfs, long nvalid,
                        T offset, T scale, int device);
  sutra_centroider_tcog(const sutra_centroider_tcog &centroider);
  ~sutra_centroider_tcog();

  string get_type();

  int set_threshold(T threshold);

  int get_cog(T *cube, T *intensities, T *centroids, int nvalid, int npix,
              int ntot);
  int get_cog(T *intensities, T *slopes, bool noise);
  int get_cog();
};

template <class T>
void get_centroids(int size, int threads, int blocks, int n, T *d_idata,
                   T *d_odata, T *ref, int *validx, int *validy, T *intensities,
                   T threshold, T scale, T offset, carma_device *device);

#endif  // _SUTRA_CENTROIDER_TCOG_H_
