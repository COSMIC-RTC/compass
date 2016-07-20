#ifndef _SUTRA_CENTROIDER_WCOG_H_
#define _SUTRA_CENTROIDER_WCOG_H_

#include <sutra_centroider.h>

class sutra_centroider_wcog: public sutra_centroider {
public:
  int npix;
  carma_obj<float> *d_weights;

public:
  sutra_centroider_wcog(carma_context *context, sutra_sensors *sensors, int nwfs, long nvalid,
      float offset, float scale, int device);
  sutra_centroider_wcog(const sutra_centroider& centroider);
  ~sutra_centroider_wcog();

  string
  get_type();

  int
  init_weights();
  int
  load_weights(float *weights, int ndim);

  int
  get_cog(carma_streams *streams, float *cube, float *subsum, float *centroids,
      int nvalid, int npix, int ntot);
  int
  get_cog(float *subsum, float *slopes, bool noise);
  int
  get_cog();
};
#endif // _SUTRA_CENTROIDER_WCOG_H_
