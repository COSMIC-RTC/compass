#ifndef _SUTRA_CENTROIDER_COG_H_
#define _SUTRA_CENTROIDER_COG_H_

#include <sutra_centroider.h>

class sutra_centroider_cog : public sutra_centroider {
 public:
 public:
  sutra_centroider_cog(carma_context *context, sutra_wfs *wfs, long nvalid,
                       float offset, float scale, int device);
  sutra_centroider_cog(const sutra_centroider_cog &centroider);
  ~sutra_centroider_cog();

  string get_type();

  int get_cog(carma_streams *streams, float *cube, float *subsum,
              float *centroids, int nvalid, int npix, int ntot);
  int get_cog(float *subsum, float *slopes, bool noise);
  int get_cog();
};

template <class T>
void get_centroids(int size, int threads, int blocks, int n, T *d_idata,
                   T *d_odata, T *alpha, int *validx, int *validy, T scale,
                   T offset, carma_device *device);

template <class T>
void get_centroids_async(int threads, int blocks, int n, carma_streams *streams,
                         T *d_idata, T *d_odata, T *alpha, T scale, T offset);

#endif  // _SUTRA_CENTROIDER_COG_H_
