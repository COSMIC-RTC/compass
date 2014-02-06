#ifndef _SUTRA_CENTROIDER_BPCOG_H_
#define _SUTRA_CENTROIDER_BPCOG_H_

#include <sutra_centroider.h>

using namespace std;

class sutra_centroider_bpcog : public sutra_centroider {
 public:
  int nmax;

 public:
  sutra_centroider_bpcog(carma_context *context, long nwfs, long nvalid, float offset, float scale, int device,int nmax);
  sutra_centroider_bpcog(const sutra_centroider_bpcog& centroider);
  ~sutra_centroider_bpcog();

  int init_nmax(int nmax); 
  int set_nmax(int nmax); 

  int init_bincube(sutra_wfs *wfs);

  int get_cog(float *cube,float *subsum,float *centroids, int nvalid, int npix, int ntot);
  int get_cog(sutra_wfs *wfs, carma_obj<float> *slopes);
  int get_cog(sutra_wfs *wfs);  

  int get_cog_async(carma_streams *streams, float *cube,float *subsum, float *centroids, int nvalid, int npix);
  int get_cog_async(sutra_wfs *wfs, carma_obj<float> *slopes);
  int get_cog_async(sutra_wfs *wfs);
};

#endif // _SUTRA_CENTROIDER_H_

