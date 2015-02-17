#ifndef _SUTRA_CENTROIDER_BPCOG_H_
#define _SUTRA_CENTROIDER_BPCOG_H_

#include <sutra_centroider.h>

using namespace std;

class sutra_centroider_bpcog: public sutra_centroider {
public:
  int nmax;
  carma_obj<float> *d_bpix;
  carma_obj<uint> *d_bpind;

public:
  sutra_centroider_bpcog(carma_context *context, sutra_sensors *sensors, int nwfs, long nvalid,
      float offset, float scale, int device, int nmax);
  sutra_centroider_bpcog(const sutra_centroider_bpcog& centroider);
  ~sutra_centroider_bpcog();

  string
  get_type();

  int
  init_nmax(int nmax);
  int
  set_nmax(int nmax);

  int
  get_cog(carma_streams *streams, float *cube, float *subsum, float *centroids,
      int nvalid, int npix, int ntot);
  int
  get_cog(float *slopes);
  int
  get_cog();
};

#endif // _SUTRA_CENTROIDER_H_
