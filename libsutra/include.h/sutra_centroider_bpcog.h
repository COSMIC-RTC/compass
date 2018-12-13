#ifndef _SUTRA_CENTROIDER_BPCOG_H_
#define _SUTRA_CENTROIDER_BPCOG_H_

#include <sutra_centroider.h>

class sutra_centroider_bpcog : public sutra_centroider {
 public:
  int nmax;
  carma_obj<float> *d_bpix;
  carma_obj<uint> *d_bpind;

 public:
  sutra_centroider_bpcog(carma_context *context, sutra_wfs *wfs, long nvalid,
                         float offset, float scale, int device, int nmax);
  sutra_centroider_bpcog(const sutra_centroider_bpcog &centroider);
  ~sutra_centroider_bpcog();

  string get_type();

  int init_nmax(int nmax);
  int set_nmax(int nmax);

  int get_cog(float *cube, float *intensities, float *centroids, int nvalid,
              int npix, int ntot);
  int get_cog(float *intensities, float *slopes, bool noise);
  int get_cog();
};

template <class T>
void subap_sortmax(int threads, int blocks, T *d_idata, T *d_odata,
                   unsigned int *values, int nmax, carma_device *device);
template <class T>
void subap_bpcentro(int threads, int blocks, int npix, T *d_idata,
                    unsigned int *values, T *d_odata, T scale, T offset);

#endif  // _SUTRA_CENTROIDER_H_
