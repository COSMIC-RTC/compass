#ifndef _SUTRA_CENTROIDER_WCOG_H_
#define _SUTRA_CENTROIDER_WCOG_H_

#include <sutra_centroider.h>

template <class Tin, class T>
class sutra_centroider_wcog : public sutra_centroider<Tin, T> {
 public:
  int npix;
  carma_obj<float> *d_weights;

 public:
  sutra_centroider_wcog(carma_context *context, sutra_wfs *wfs, long nvalid,
                        float offset, float scale, bool filter_TT, int device);
  sutra_centroider_wcog(const sutra_centroider_wcog &centroider);
  ~sutra_centroider_wcog();

  string get_type();

  int set_npix(int npix);
  int init_weights();
  int load_weights(float *weights, int ndim);

  int get_cog(float *cube, float *intensities, T *centroids, int nvalid, int npix,
              int ntot);
  int get_cog(float *intensities, T *slopes, bool noise);
  int get_cog();
};

template <class T>
void get_centroids(int size, int threads, int blocks, int n, float *d_idata,
                   T *d_odata, T *ref, int *validx, int *validy, float *intensities,
                   float *d_weights, float scale, float offset,
                   carma_device *device);

template <class T>
int fillweights(T *d_out, T *d_in, int npix, int N, carma_device *device);
#endif  // _SUTRA_CENTROIDER_WCOG_H_
