#ifndef _SUTRA_CENTROIDER_COG_H_
#define _SUTRA_CENTROIDER_COG_H_

#include <sutra_centroider.h>

template <class Tin, class T>
class sutra_centroider_cog : public sutra_centroider<Tin, T> {
 public:
  sutra_centroider_cog(carma_context *context, sutra_wfs *wfs, long nvalid,
                       float offset, float scale, int device);
  sutra_centroider_cog(const sutra_centroider_cog &centroider);
  ~sutra_centroider_cog();

  string get_type();

  int get_cog(float *cube, float *intensities, T *centroids, int nvalid,
              int npix, int ntot);
  int get_cog(float *intensities, T *slopes, bool noise);
  int get_cog();
};

template <class T>
void get_centroids(int size, int threads, int blocks, int n, float *d_idata,
                   T *d_odata, T *ref, int *validx, int *validy,
                   float *intensities, float scale, float offset,
                   carma_device *device);

template <class T>
void get_centroids_async(int threads, int blocks, int n, carma_streams *streams,
                         T *d_idata, T *d_odata, T *alpha, float scale,
                         float offset);

#endif  // _SUTRA_CENTROIDER_COG_H_
