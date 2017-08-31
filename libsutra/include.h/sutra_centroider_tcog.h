#ifndef _SUTRA_CENTROIDER_TCOG_H_
#define _SUTRA_CENTROIDER_TCOG_H_

#include <sutra_centroider.h>

class sutra_centroider_tcog: public sutra_centroider {
 public:
  float threshold;

 public:
  sutra_centroider_tcog(carma_context *context, sutra_sensors *sensors, int nwfs, long nvalid,
                        float offset, float scale, int device);
  sutra_centroider_tcog(const sutra_centroider& centroider);
  ~sutra_centroider_tcog();

  string
  get_type();

  int
  set_threshold(float threshold);

  int
  get_cog(carma_streams *streams, float *cube, float *subsum, float *centroids,
          int nvalid, int npix, int ntot);
  int
  get_cog(float *subsum, float *slopes, bool noise);
  int
  get_cog();
};

#endif // _SUTRA_CENTROIDER_TCOG_H_
