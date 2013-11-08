#ifndef _SUTRA_CENTROIDER_H_
#define _SUTRA_CENTROIDER_H_

#include <sutra_wfs.h>

using namespace std;

class sutra_centroider {
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

  carma_obj<float>           *d_bincube;
  carma_obj<float>           *d_subsum;

  carma_obj<float>          *d_centroids;     
  carma_obj<float>          *d_refcentro;

  // this is for centroiding
  carma_obj<float>          *d_validpix;     
  carma_obj<int>            *d_validindx;

  carma_obj<float>          *d_weights; 
    
  carma_obj<cuFloatComplex> *d_corrfnct;
  carma_obj<cuFloatComplex> *d_corrspot;
  carma_obj<float>          *d_corrnorm;     
  carma_obj<int>            *d_corrmax;     
  carma_obj<float>          *d_corr;     
  carma_obj<float>          *d_interpmat;     

  carma_context *current_context;

 public:
  sutra_centroider(carma_context *context, long nwfs, long nvalid, float offset, float scale, int device, const char* typec);
  sutra_centroider(const sutra_centroider& centroider);
  ~sutra_centroider();

  int init_nmax(int nmax); 
  int set_nmax(int nmax); 
  int set_threshold(float threshold); 

  int init_bincube(sutra_wfs *wfs);

  int init_weights(sutra_wfs *wfs);
  int load_weights(float *weights, int ndim);

  int init_corr(sutra_wfs *wfs, int isizex, int isizey, float *interpmat);
  int load_corr(float *corr, float *corr_norm, int ndim);
 
  int get_cog(float *cube,float *subsum,float *centroids, int nvalid, int npix, int ntot);
  int get_cog(sutra_wfs *wfs, carma_obj<float> *slopes);
  int get_cog(sutra_wfs *wfs);  

  int get_cog_async(carma_streams *streams, float *cube,float *subsum, float *centroids, int nvalid, int npix);
  int get_cog_async(sutra_wfs *wfs, carma_obj<float> *slopes);
  int get_cog_async(sutra_wfs *wfs);

  int get_tcog(float threshold,float *cube,float *subsum, float *centroids, int nvalid, int npix, int ntot);
  int get_tcog(float threshold,sutra_wfs *wfs, float *slopes);
  int get_tcog(float threshold,sutra_wfs *wfs);

  int get_bpcog(int nmax, int npix, int nvalid, float *cube, float *centroids);
  int get_bpcog(int nmax,sutra_wfs *wfs, float *slopes);
  int get_bpcog(int nmax,sutra_wfs *wfs);

  int get_wcog(float *weights,float *cube,float *subsum,float *centroids, int nvalid, int npix, int ntot);
  int get_wcog(sutra_wfs *wfs, float *slopes);
  int get_wcog(sutra_wfs *wfs);

  int get_nmax(float *cube, int npix, int nvalid,int ntot);
  int get_nmax(sutra_wfs *wfs);

  int get_corr(sutra_wfs *wfs, float *slopes);
  int get_corr(sutra_wfs *wfs);

  int get_pyr(float *cube,float *subsum, float *centroids, int *subindx, int *subindy, int nvalid, int ns, int nim);
  int get_pyr(sutra_wfs *wfs, carma_obj<float> *slopes);
  int get_pyr(sutra_wfs *wfs);  
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
template <class T> void get_centroids_async(int threads, int blocks, int n, carma_streams *streams, T *d_idata, T *d_odata,  T *alpha, T scale, T offset);
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
#endif // _SUTRA_CENTROIDER_H_

