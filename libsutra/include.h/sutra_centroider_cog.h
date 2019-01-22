#ifndef _SUTRA_CENTROIDER_COG_H_
#define _SUTRA_CENTROIDER_COG_H_

#include <sutra_centroider.h>

template <typename T>
class sutra_centroider_cog : public sutra_centroider<T> {
 public:
 public:
  sutra_centroider_cog(carma_context *context, sutra_wfs *wfs, long nvalid,
                       T offset, T scale, int device);
  sutra_centroider_cog(const sutra_centroider_cog &centroider);
  ~sutra_centroider_cog();

  string get_type();

  int get_cog(T *cube, T *intensities, T *centroids, int nvalid, int npix,
              int ntot);
  int get_cog(T *intensities, T *slopes, bool noise);
  int get_cog();
};

template <class T>
void get_centroids(int size, int threads, int blocks, int n, T *d_idata,
                   T *d_odata, T *ref, int *validx, int *validy, T *intensities,
                   T scale, T offset, carma_device *device);

template <class T>
void get_centroids_async(int threads, int blocks, int n, carma_streams *streams,
                         T *d_idata, T *d_odata, T *alpha, T scale, T offset);

#endif  // _SUTRA_CENTROIDER_COG_H_
