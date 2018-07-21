#ifndef _SUTRA_CENTROIDER_CORR_H_
#define _SUTRA_CENTROIDER_CORR_H_

#include <sutra_centroider.h>

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
  sutra_centroider_corr(carma_context *context, sutra_wfs *wfs, long nvalid,
                        float offset, float scale, int device);
  sutra_centroider_corr(const sutra_centroider_corr &centroider);
  ~sutra_centroider_corr();

  string get_type();

  int init_corr(int isizex, int isizey, float *interpmat);
  int load_corr(float *corr, float *corr_norm, int ndim);

  int init_bincube(int npix);

  int get_cog(carma_streams *streams, float *cube, float *subsum,
              float *centroids, int nvalid, int npix, int ntot);
  int get_cog(float *subsum, float *slopes, bool noise);
  int get_cog();
};

#endif  // _SUTRA_CENTROIDER_CORR_H_
