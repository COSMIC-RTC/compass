#ifndef _SUTRA_CENTROIDER_H_
#define _SUTRA_CENTROIDER_H_

#include <sutra_wfs.h>
#include <string>

using namespace std;

class sutra_centroider {
public:
  int device;
  int nwfs;
  int nvalid;

  float offset;
  float scale;

  carma_context *current_context;

public:
  virtual ~sutra_centroider()=0;

  bool is_type(string typec){
    return (typec.compare(get_type())==0);
  }

  virtual string get_type()=0;

  virtual int init_bincube(sutra_wfs *wfs)=0;

  virtual int get_cog(carma_streams *streams, float *cube, float *subsum, float *centroids, int nvalid,
      int npix, int ntot)=0;
  virtual int get_cog(sutra_wfs *wfs, float *slopes)=0;
  virtual int get_cog(sutra_wfs *wfs)=0;
};

int fillweights(float *d_out, float *d_in, int npix, int N, int device);
int fillcorr(cuFloatComplex *d_out, float *d_in, int npix_in, int npix_out,
    int N, int nvalid, int device);
int correl(cuFloatComplex *d_odata, cuFloatComplex *d_idata, int N, int device);
int roll2real(float *d_odata, cuFloatComplex *d_idata, int n, int Npix, int N,
    int device);
int corr_norm(float *d_odata, float *d_idata, int Npix, int N, int device);
int convert_centro(float *d_odata, float *d_idata, float offset, float scale,
    int N, int device);
//int fillval_corr(cuFloatComplex *d_out, float val, int npix_in, int npix_out, int N, int device);

// CUDA templates
template<class T> void get_centroids(int size, int threads, int blocks, int n,
    T *d_idata, T *d_odata, T *alpha, T scale, T offset);
template<class T> void get_centroids_async(int threads, int blocks, int n,
    carma_streams *streams, T *d_idata, T *d_odata, T *alpha, T scale,
    T offset);
template<class T> void get_centroids(int size, int threads, int blocks, int n,
    T *d_idata, T *d_odata, T *alpha, T thresh, T scale, T offset);
template<class T> void get_centroids(int size, int threads, int blocks, int n,
    T *d_idata, T *d_odata, T *alpha, T *weights, T scale, T offset);
template<class T> void subap_centromax(int threads, int blocks, T *d_idata,
    T *d_odata, int npix, int nmax, T scale, T offset);
template<class T> void subap_sortmax(int size, int threads, int blocks,
    T *d_idata, T *d_odata, int *values, int nmax);
template<class T> void subap_sortmaxi(int threads, int blocks, T *d_idata,
    int *values, int nmax, int offx, int offy, int npix, int Npix);
template<class T> void subap_pinterp(int threads, int blocks, T *d_idata,
    int *values, T *d_centroids, T *d_matinterp, int sizex, int sizey,
    int nvalid, int Npix, T scale, T offset);
template<class T> void pyr_slopes(T *d_odata, T *d_idata, int *subindx,
    int *subindy, T *subsum, int ns, int nvalid, int nim, int device);
#endif // _SUTRA_CENTROIDER_H_
