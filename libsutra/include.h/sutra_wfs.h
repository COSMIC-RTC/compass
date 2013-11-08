#ifndef _SUTRA_WFS_H_
#define _SUTRA_WFS_H_

#include <vector>
#include <map>
#include <sutra_telemetry.h>
#include <sutra_target.h>
#include <sutra_phase.h>
#include <sutra_lgs.h>
//#include <sutra_slopes.h>

extern "C" {
    void my_abort(int err);
}

using namespace std;

class sutra_wfs {
 public:
  int                      device;
  string                   type;
  long                     nxsub;
  long                     nvalid;
  long                     npix;
  long                     nrebin;
  long                     nfft;
  long                     ntot;
  long                     npup;
  long                     nphase;
  long                     nmaxhr;
  long                     nffthr;
  float                    subapd;
  float                    nphot;
  float                    noise;
  bool                     lgs;
  bool                     kernconv;

  carma_obj<cuFloatComplex> *d_camplipup;
  carma_obj<cuFloatComplex> *d_camplifoc;
  carma_obj<cuFloatComplex> *d_fttotim;
  carma_obj<cuFloatComplex> *d_ftkernel;

  carma_obj<float>          *d_pupil;
  carma_obj<float>          *d_hrimg;
  carma_obj<float>          *d_bincube;
  carma_obj<float>          *d_binimg;
  carma_obj<float>          *d_subsum;
  carma_obj<float>          *d_offsets;
  carma_obj<float>          *d_fluxPerSub;
  carma_obj<float>          *d_sincar;
  carma_obj<float>          *d_submask;
  carma_obj<int>            *d_hrmap;

  carma_obj<int>            *d_isvalid;    // nxsub x nxsub
  carma_obj<float>          *d_slopes;

  carma_host_obj<float>     *image_telemetry;

  // sh only
  carma_obj<int>            *d_phasemap;
  carma_obj<int>            *d_binmap;
  carma_obj<int>            *d_validsubsx;  // nvalid
  carma_obj<int>            *d_validsubsy;  // nvalid
  carma_obj<int>            *d_istart;     // nxsub 
  carma_obj<int>            *d_jstart;     // nxsub

  // pyramid only
  carma_obj<float>          *d_psum;
  carma_obj<cuFloatComplex> *d_phalfxy;
  carma_obj<cuFloatComplex> *d_poffsets;

  carma_host_obj<int>       *pyr_cx;
  carma_host_obj<int>       *pyr_cy;

  sutra_source              *d_gs;

  carma_streams 		   *streams;
  int  			   nstreams;

  carma_context *current_context;

 public:
  sutra_wfs(carma_context *context, const char* type, long nxsub, long nvalid, long npix, long nphase, 
	   long nrebin, long nfft, long ntot, long npup,float pdiam,float nphotons, int lgs, int device);
  sutra_wfs(const sutra_wfs& wfs);
  ~sutra_wfs();

  int wfs_initarrays(int *phasemap,int *hrmap, int *binmap,float *offsets, 
			     float *pupil, float *fluxPerSub, int *isvalid, int *validsubsx, int *validsubsy, 
			     int *istart, int *jstart, cuFloatComplex *kernel);
  int wfs_initarrays(cuFloatComplex *halfxy,cuFloatComplex *offsets, float *focmask, 
		     float *pupil, int *isvalid, int *cx, int *cy, float *sincar,int *phasemap,
		     int *validsubsx, int *validsubsy);
  int wfs_initgs(float xpos,float ypos,float lambda, float mag, long size,float noise, long seed);
  int load_kernels(float *lgskern);
  int sensor_trace(sutra_atmos *yatmos);
  int sensor_trace(sutra_dms *ydm, int rst);
  int comp_image_tele();
  int comp_image();
  int slopes_geom(int type, float *slopes);
  int slopes_geom(int type);

private:
  int comp_sh_generic();
  int comp_pyr_generic();
};

class sutra_sensors {
 public:
  int                 nsensors;
  vector<sutra_wfs *>  d_wfs;
     
 public:
  sutra_sensors(carma_context *context, const char* type, int nwfs,long *nxsub,long *nvalid,long *npix,
	       long *nphase, long *nrebin,long *nfft, long *ntot, long npup, float *pdiam, float *nphot,
	       int *lgs, int device);
  ~sutra_sensors();

  int sensors_initgs(float *xpos,float *ypos,float *lambda, float *mag, long *size, float *noise, long *seed);
  int sensors_initgs(float *xpos,float *ypos,float *lambda, float *mag, long *size, float *noise);
  int sensors_initgs(float *xpos,float *ypos,float *lambda, float *mag, long *size);
};

// General utilities
int fillcamplipup(cuFloatComplex *amplipup, float *phase, float *offset, float *mask, float scale, int *istart, int *jstart, 
		   int *ivalid, int *jvalid,int nphase, int npup, int Nfft, int Ntot, int device);
int indexfill(cuFloatComplex *d_odata,cuFloatComplex *d_idata,int *indx,int nx, int Nx, int N,int device);
int fillbincube(float *bcube,cuFloatComplex  *hrimage, int *indxpix, int Nfft, int Npix, int Nrebin, int Nsub, 
		int device);
int fillbincube_async(carma_streams *streams, float *bcube,cuFloatComplex  *hrimage, int *indxpix, int Nfft, 
		      int Npix, int Nrebin, int Nsub, int device);
int fillbinimg(float *bimage, float *bcube, int npix, int nsub, int Nsub, int *ivalid, int *jvalid, bool add, 
	       int device);
int fillbinimg_async(carma_streams *streams, carma_obj<float> *bimage, carma_obj<float> *bcube, int npix, int nsub, 
		     int Nsub, int *ivalid, int *jvalid, bool add, int device);
int fillbinimg_async(carma_host_obj<float> *image_telemetry, float *bimage, float *bcube, int npix,
		     int nsub, int Nsub, int *ivalid, int *jvalid,  int nim, bool add, int device);
int convolve(cuFloatComplex *d_odata,cuFloatComplex *d_idata,int N,int device);
int convolve_cube(cuFloatComplex *d_odata,cuFloatComplex *d_idata,int N,int n, int device);

// CUDA templates
// this is for cog
template <class T> void subap_reduce(int size, int threads, int blocks, T *d_idata, T *d_odata);
template <class T> void subap_reduce_async(int threads, int blocks, carma_streams *streams, T *d_idata, T *d_odata);
// this is for tcog
template <class T> void subap_reduce(int size, int threads, int blocks, T *d_idata, T *d_odata, T thresh);
// this is for wcog
template <class T> void subap_reduce(int size, int threads, int blocks, T *d_idata, T *d_odata, T *weights);
template <class T> void phase_reduce(int threads, int blocks, T *d_idata, T *d_odata, int *indx, T alpha);
template <class T> void phase_derive(int size, int threads, int blocks, int n, T *d_idata, T *d_odata, 
				     int *indx, T *mask, T alpha, float *fluxPerSub);
template <class Tout, class Tin> void pyr_getpup(Tout *d_odata, Tin *d_idata, Tout *d_offsets, Tin *d_pup, int np, 
						 int device);
template <class T> void pyr_rollmod(T *d_odata, T *d_idata, T *d_mask, float cx, float cy, int np, int ns, int device);
template <class T> void pyr_fillbin(T *d_odata, T *d_idata, int nrebin, int np, int ns, int nim, int device);
template <class Tin, class Tout> void pyr_abs2(Tout *d_odata, Tin *d_idata, Tout fact, int ns, int nim, int device);
template <class Tout, class Tin> void pyr_submask(Tout *d_odata, Tin *d_mask, int n, int device);
template <class Tout,class Tin> void pyr_abs(Tout *d_odata, Tin *d_idata, int ns, int nim, int device);
template <class Tout, class Tin> void pyr_submask3d(Tout *d_odata, Tin *d_mask, int n, int nim, int device);
template <class T> void pyr_subsum(T *d_odata, T *d_idata, int *subindx, int *subindy, int ns, int nvalid, int nim, 
				   int device);
template <class T> void pyr_fact(T *d_data, T fact, int n, int nim, int device);
void pyr_fact(cuFloatComplex *d_data, float fact, int n, int nim, int device);
void pyr_fact(float *d_data, float fact1, float *fact2, int n, int nim, int device);
#endif // _SUTRA_WFS_H_
