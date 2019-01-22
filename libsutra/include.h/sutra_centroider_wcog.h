#ifndef _SUTRA_CENTROIDER_WCOG_H_
#define _SUTRA_CENTROIDER_WCOG_H_

#include <sutra_centroider.h>

template <typename T>
class sutra_centroider_wcog : public sutra_centroider<T> {
 public:
  int npix;
  carma_obj<T> *d_weights;

 public:
  sutra_centroider_wcog(carma_context *context, sutra_wfs *wfs, long nvalid,
                        T offset, T scale, int device);
  sutra_centroider_wcog(const sutra_centroider_wcog &centroider);
  ~sutra_centroider_wcog();

  string get_type();

  int set_npix(int npix);
  int init_weights();
  int load_weights(T *weights, int ndim);

  int get_cog(T *cube, T *intensities, T *centroids, int nvalid, int npix,
              int ntot);
  int get_cog(T *intensities, T *slopes, bool noise);
  int get_cog();
};

template <class T>
void get_centroids(int size, int threads, int blocks, int n, T *d_idata,
                   T *d_odata, T *ref, int *validx, int *validy, T *intensities,
                   T *d_weights, T scale, T offset, carma_device *device);

template <class T>
int fillweights(T *d_out, T *d_in, int npix, int N, carma_device *device);
#endif  // _SUTRA_CENTROIDER_WCOG_H_
