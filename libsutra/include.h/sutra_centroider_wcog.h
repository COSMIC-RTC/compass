#ifndef _SUTRA_CENTROIDER_WCOG_H_
#define _SUTRA_CENTROIDER_WCOG_H_

#include <sutra_centroider.h>

using namespace std;

class sutra_centroider_wcog : public sutra_centroider {
 public:
  int npix;
  carma_obj<float> *d_weights;

 public:
  sutra_centroider_wcog(carma_context *context, long nwfs, long nvalid, float offset, float scale, int device);
  sutra_centroider_wcog(const sutra_centroider& centroider);
  ~sutra_centroider_wcog();

  int init_bincube(sutra_wfs *wfs);

  int init_weights(sutra_wfs *wfs);
  int load_weights(float *weights, int ndim);

  int get_cog(float *cube,float *subsum,float *centroids, int nvalid, int npix, int ntot);
  int get_cog(sutra_wfs *wfs, carma_obj<float> *slopes);
  int get_cog(sutra_wfs *wfs);  

  int get_cog_async(carma_streams *streams, float *cube,float *subsum, float *centroids, int nvalid, int npix);
  int get_cog_async(sutra_wfs *wfs, carma_obj<float> *slopes);
  int get_cog_async(sutra_wfs *wfs);

};
#endif // _SUTRA_CENTROIDER_WCOG_H_

