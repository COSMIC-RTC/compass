#ifndef _SUTRA_CENTROIDER_H_
#define _SUTRA_CENTROIDER_H_

#include <sutra_acquisim.h>
#include <sutra_wfs.h>
#include <string>

class sutra_centroider {
 public:
  int device;
  sutra_wfs *wfs;
  int nvalid;

  float offset;
  float scale;

  carma_context *current_context;

  carma_obj<float> *d_bincube;
  carma_obj<float> *d_img;
  carma_obj<int> *d_validx;
  carma_obj<int> *d_validy;

 protected:
  sutra_centroider(carma_context *context, sutra_wfs *wfs, long nvalid,
                   float offset, float scale, int device);

 public:
  virtual ~sutra_centroider();
  int load_validpos(int *ivalid, int *jvalid, int N);
  int fill_bincube(int npix);
  int load_img(float *img, int n);
  int load_pyrimg(float *img, int n);
  bool is_type(string typec) { return (typec.compare(get_type()) == 0); }
  int load_img_gpu(float *img, int n);

  virtual string get_type() = 0;

  virtual int get_cog(carma_streams *streams, float *cube, float *subsum,
                      float *centroids, int nvalid, int npix, int ntot) = 0;
  virtual int get_cog(float *subsum, float *slopes, bool noise) = 0;
  virtual int get_cog() = 0;
};

int fillweights(float *d_out, float *d_in, int npix, int N,
                carma_device *device);
int fillcorr(cuFloatComplex *d_out, float *d_in, int npix_in, int npix_out,
             int N, int nvalid, carma_device *device);
int correl(cuFloatComplex *d_odata, cuFloatComplex *d_idata, int N,
           carma_device *device);
int roll2real(float *d_odata, cuFloatComplex *d_idata, int n, int Npix, int N,
              carma_device *device);
int corr_norm(float *d_odata, float *d_idata, int Npix, int N,
              carma_device *device);
int convert_centro(float *d_odata, float *d_idata, float offset, float scale,
                   int N, carma_device *device);
// int fillval_corr(cuFloatComplex *d_out, float val, int npix_in, int npix_out,
// int N, carma_device *device);

// CUDA templates
template <class T>
void get_centroids(int size, int threads, int blocks, int n, T *d_idata,
                   T *d_odata, T *alpha, T scale, T offset,
                   carma_device *device);
template <class T>
void get_centroids_async(int threads, int blocks, int n, carma_streams *streams,
                         T *d_idata, T *d_odata, T *alpha, T scale, T offset);
template <class T>
void get_centroids(int size, int threads, int blocks, int n, T *d_idata,
                   T *d_odata, T *alpha, T thresh, T scale, T offset,
                   carma_device *device);
template <class T>
void get_centroids(int size, int threads, int blocks, int n, T *d_idata,
                   T *d_odata, T *alpha, T *weights, T scale, T offset,
                   carma_device *device);
template <class T>
void subap_bpcentro(int threads, int blocks, int npix, T *d_idata,
                    unsigned int *values, T *d_odata, T scale, T offset);
template <class T>
void subap_centromax(int threads, int blocks, T *d_idata, T *d_odata, int npix,
                     int nmax, T scale, T offset);
template <class T>
void subap_centromax2(int threads, int blocks, T *d_idata, T *d_odata,
                      T *d_minim, int npix, int nmax, T scale, T offset);
template <class T>
void subap_sortmax(int threads, int blocks, T *d_idata, T *d_odata,
                   unsigned int *values, int nmax, carma_device *device);
template <class T>
void subap_sortmaxi(int threads, int blocks, T *d_idata, int *values, int nmax,
                    int offx, int offy, int npix, int Npix);
template <class T>
void subap_pinterp(int threads, int blocks, T *d_idata, int *values,
                   T *d_centroids, T *d_matinterp, int sizex, int sizey,
                   int nvalid, int Npix, T scale, T offset);
template <class T>
void pyr_slopes(T *d_odata, T *d_idata, int *subindx, int *subindy, T *subsum,
                int ns, int nvalid, int nim, carma_device *device);

template <class T>
void pyr2_slopes(T *d_odata, T *d_idata, int *subindx, int *subindy, T *subsum,
                 int ns, int nvalid, T scale, T valid_thresh, int do_sin,
                 carma_device *device);

template <class T>
void roof_slopes(T *d_odata, T *d_idata, int *subindx, int *subindy, T *subsum,
                 int ns, int nvalid, int nim, carma_device *device);

#endif  // _SUTRA_CENTROIDER_H_
