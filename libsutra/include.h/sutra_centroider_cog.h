#ifndef _SUTRA_CENTROIDER_COG_H_
#define _SUTRA_CENTROIDER_COG_H_

#include <sutra_centroider.h>

class sutra_centroider_cog : public sutra_centroider {
 public:
 public:
  sutra_centroider_cog(carma_context *context, sutra_sensors *sensors, int nwfs,
                       long nvalid, float offset, float scale, int device);
  sutra_centroider_cog(const sutra_centroider_cog &centroider);
  ~sutra_centroider_cog();

  string get_type();

  int get_cog(carma_streams *streams, float *cube, float *subsum,
              float *centroids, int nvalid, int npix, int ntot);
  int get_cog(float *subsum, float *slopes, bool noise);
  int get_cog();
};

#endif  // _SUTRA_CENTROIDER_COG_H_
