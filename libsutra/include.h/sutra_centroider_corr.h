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

  int set_npix(int npix);

  int get_cog(float *cube, float *intensities, float *centroids, int nvalid,
              int npix, int ntot);
  int get_cog(float *intensities, float *slopes, bool noise);
  int get_cog();
};

template <class T>
void subap_sortmaxi(int threads, int blocks, T *d_idata, int *values, int nmax,
                    int offx, int offy, int npix, int Npix);
template <class T>
void subap_pinterp(int threads, int blocks, T *d_idata, int *values,
                   T *d_centroids, T *d_matinterp, int sizex, int sizey,
                   int nvalid, int Npix, T scale, T offset);

template <class Tcu, class T>
int fillcorr(Tcu *d_out, T *d_in, int npix_in, int npix_out, int N, int nvalid,
             carma_device *device);

template <class T>
int correl(T *d_odata, T *d_idata, int N, carma_device *device);

template <class Tcu, class T>
int roll2real(T *d_odata, Tcu *d_idata, int n, int Npix, int N,
              carma_device *device);

template <class T>
int corr_norm(T *d_odata, T *d_idata, int Npix, int N, carma_device *device);

#endif  // _SUTRA_CENTROIDER_CORR_H_
