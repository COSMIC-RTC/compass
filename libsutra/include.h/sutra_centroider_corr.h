#ifndef _SUTRA_CENTROIDER_CORR_H_
#define _SUTRA_CENTROIDER_CORR_H_

#include <sutra_centroider.h>

using namespace std;

class sutra_centroider_corr : public sutra_centroider {
 public:
  int npix;
  int interp_sizex;
  int interp_sizey;
  carma_obj<cuFloatComplex> *d_corrfnct;
  carma_obj<cuFloatComplex> *d_corrspot;
  carma_obj<float> *d_corrnorm;
  carma_obj<int> *d_corrmax;
  carma_obj<float> *d_corr;
  carma_obj<float> *d_interpmat;

 public:
  sutra_centroider_corr(carma_context *context, long nwfs, long nvalid, float offset, float scale, int device);
  sutra_centroider_corr(const sutra_centroider_corr& centroider);
  ~sutra_centroider_corr();

  string get_type();

  int init_corr(sutra_wfs *wfs, int isizex, int isizey, float *interpmat);
  int load_corr(float *corr, float *corr_norm, int ndim);

  int init_bincube(sutra_wfs *wfs);

  int get_cog(carma_streams *streams, float *cube, float *subsum, float *centroids, int nvalid,
      int npix, int ntot);
  int get_cog(sutra_wfs *wfs, float *slopes);
  int get_cog(sutra_wfs *wfs);
};

#endif // _SUTRA_CENTROIDER_CORR_H_

