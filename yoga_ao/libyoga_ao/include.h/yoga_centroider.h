#ifndef _YOGA_CENTROIDER_H_
#define _YOGA_CENTROIDER_H_

#include <yoga_wfs.h>

using namespace std;

class yoga_centroider {
 public:
  string                   typec;
  int                      device;
  int                      nwfs;
  int                      nvalid;
  int                      npix;
  float                    pixsize;
  int                      interp_sizex;
  int                      interp_sizey;

  int                      nmax;
  float                    threshold;

  float                    offset;
  float                    scale;

  yoga_obj<float>           *d_bincube;
  yoga_obj<float>           *d_subsum;

  yoga_obj<float>          *d_centroids;     
  yoga_obj<float>          *d_refcentro;

  // this is for centroiding
  yoga_obj<float>          *d_validpix;     
  yoga_obj<int>            *d_validindx;

  yoga_obj<float>          *d_weights; 
    
  yoga_obj<cuFloatComplex> *d_corrfnct;
  yoga_obj<cuFloatComplex> *d_corrspot;
  yoga_obj<float>          *d_corrnorm;     
  yoga_obj<int>            *d_corrmax;     
  yoga_obj<float>          *d_corr;     
  yoga_obj<float>          *d_interpmat;     

  yoga_context *current_context;

 public:
  yoga_centroider(yoga_context *context, long nwfs, long nvalid, float offset, float scale, int device, const char* typec);
  yoga_centroider(const yoga_centroider& centroider);
  ~yoga_centroider();

  int init_nmax(int nmax); 
  int set_nmax(int nmax); 
  int set_threshold(float threshold); 

  int init_bincube(yoga_wfs *wfs);

  int init_weights(yoga_wfs *wfs);
  int load_weights(float *weights, int ndim);

  int init_corr(yoga_wfs *wfs, int isizex, int isizey, float *interpmat);
  int load_corr(float *corr, float *corr_norm, int ndim);
 
  int get_cog(float *cube,float *subsum,float *centroids, int nvalid, int npix, int ntot);
  int get_cog(yoga_wfs *wfs, yoga_obj<float> *slopes);
  int get_cog(yoga_wfs *wfs);  

  int get_cog_async(yoga_streams *streams, float *cube,float *subsum, float *centroids, int nvalid, int npix);
  int get_cog_async(yoga_wfs *wfs, yoga_obj<float> *slopes);
  int get_cog_async(yoga_wfs *wfs);

  int get_tcog(float threshold,float *cube,float *subsum, float *centroids, int nvalid, int npix, int ntot);
  int get_tcog(float threshold,yoga_wfs *wfs, float *slopes);
  int get_tcog(float threshold,yoga_wfs *wfs);

  int get_bpcog(int nmax, int npix, int nvalid, float *cube, float *centroids);
  int get_bpcog(int nmax,yoga_wfs *wfs, float *slopes);
  int get_bpcog(int nmax,yoga_wfs *wfs);

  int get_wcog(float *weights,float *cube,float *subsum,float *centroids, int nvalid, int npix, int ntot);
  int get_wcog(yoga_wfs *wfs, float *slopes);
  int get_wcog(yoga_wfs *wfs);

  int get_nmax(float *cube, int npix, int nvalid,int ntot);
  int get_nmax(yoga_wfs *wfs);

  int get_corr(yoga_wfs *wfs, float *slopes);
  int get_corr(yoga_wfs *wfs);

  int get_pyr(float *cube,float *subsum, float *centroids, int *subindx, int *subindy, int nvalid, int ns, int nim);
  int get_pyr(yoga_wfs *wfs, yoga_obj<float> *slopes);
  int get_pyr(yoga_wfs *wfs);  
};

int fillweights(float *d_out, float *d_in, int npix, int N, int device);
int fillcorr(cuFloatComplex *d_out, float *d_in, int npix_in, int npix_out, int N, int nvalid, int device);
int correl(cuFloatComplex *d_odata,cuFloatComplex *d_idata,int N,int device);
int roll2real(float *d_odata,cuFloatComplex *d_idata, int n, int Npix, int N, int device);
int corr_norm(float *d_odata,float *d_idata,int Npix, int N,int device);
int convert_centro(float *d_odata,float *d_idata,float offset, float scale, int N,int device);
//int fillval_corr(cuFloatComplex *d_out, float val, int npix_in, int npix_out, int N, int device);

// CUDA templates
template <class T> void get_centroids(int size, int threads, int blocks, int n, T *d_idata, T *d_odata, 
				      T *alpha, T scale, T offset);
template <class T> void get_centroids_async(int threads, int blocks, int n, yoga_streams *streams, T *d_idata, T *d_odata,  T *alpha, T scale, T offset);
template <class T> void get_centroids(int size, int threads, int blocks, int n, T *d_idata, T *d_odata, 
				      T *alpha, T thresh, T scale, T offset);
template <class T> void get_centroids(int size, int threads, int blocks, int n, T *d_idata, T *d_odata, 
				      T *alpha, T *weights, T scale, T offset);
template <class T> void subap_centromax(int threads, int blocks, T *d_idata, T *d_odata, int npix, int nmax
					, T scale, T offset);
template <class T> void subap_sortmax(int size, int threads, int blocks, T *d_idata, T *d_odata, int *values, 
				      int nmax);
template <class T> void subap_sortmaxi(int threads, int blocks, T *d_idata,  int *values, int nmax,
				       int offx, int offy, int npix, int Npix);
template <class T> void subap_pinterp(int threads, int blocks, T *d_idata,  int *values, T *d_centroids,
				      T *d_matinterp, int sizex,int sizey, int nvalid, int Npix, T scale, T offset);
template <class T> void pyr_slopes(T *d_odata, T *d_idata, int *subindx, int *subindy, T *subsum, int ns, int nvalid, 
				   int nim, int device);
#endif // _YOGA_CENTROIDER_H_

