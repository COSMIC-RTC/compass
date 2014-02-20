#ifndef _SUTRA_CENTROIDER_PYR_H_
#define _SUTRA_CENTROIDER_PYR_H_

#include <sutra_centroider.h>

using namespace std;

class sutra_centroider_pyr : public sutra_centroider {
 public:

 public:
  sutra_centroider_pyr(carma_context *context, long nwfs, long nvalid, float offset, float scale, int device);
  sutra_centroider_pyr(const sutra_centroider_pyr& centroider);
  ~sutra_centroider_pyr();

  string get_type();

  int get_pyr(float *cube,float *subsum, float *centroids, int *subindx, int *subindy, int nvalid, int ns, int nim);

  int init_bincube(sutra_wfs *wfs);

  int get_cog(carma_streams *streams, float *cube, float *subsum, float *centroids, int nvalid,
      int npix, int ntot);
  int get_cog(sutra_wfs *wfs, float *slopes);
  int get_cog(sutra_wfs *wfs);
};

#endif // _SUTRA_CENTROIDER_PYR_H_

