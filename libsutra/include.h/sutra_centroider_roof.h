#ifndef _SUTRA_CENTROIDER_ROOF_H_
#define _SUTRA_CENTROIDER_ROOF_H_

#include <sutra_centroider.h>

class sutra_centroider_roof: public sutra_centroider {
public:

public:
  sutra_centroider_roof(carma_context *context, sutra_sensors *sensors, int nwfs, long nvalid,
      float offset, float scale, int device);
  sutra_centroider_roof(const sutra_centroider_roof& centroider);
  ~sutra_centroider_roof();

  string
  get_type();

  int
  get_roof(float *cube, float *subsum, float *centroids, int *subindx,
      int *subindy, int nvalid, int ns, int nim);

  int
  get_cog(carma_streams *streams, float *cube, float *subsum, float *centroids,
      int nvalid, int npix, int ntot);
  int
  get_cog(float *subsum, float *slopes, bool noise);
  int
  get_cog();
};

#endif // _SUTRA_CENTROIDER_ROOF_H_
