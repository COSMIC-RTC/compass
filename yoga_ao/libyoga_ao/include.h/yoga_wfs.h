#ifndef _YOGA_WFS_H_
#define _YOGA_WFS_H_

#include <vector>
#include <map>
#include <yoga_telemetry.h>
#include <yoga_target.h>
#include <yoga_phase.h>
#include <yoga_lgs.h>
//#include <yoga_slopes.h>

extern "C" {
    void my_abort(int err);
}

using namespace std;

class yoga_wfs_child {
 public:
  int                      device;
  string                   type;
  long                     nvalid;
  long                     nfft;
  long                     ntot;
  long                     npup;

  yoga_obj<cuFloatComplex> *d_camplipup;
  yoga_obj<cuFloatComplex> *d_camplifoc;
  yoga_obj<cuFloatComplex> *d_fttotim;
  yoga_obj<float>          *d_hrimg;
  yoga_obj<float>          *d_bincube;
  yoga_obj<float>          *d_submask;

 public:
  yoga_wfs_child(const char* type, long nvalid, long nfft, long ntot, long npup, int device);
  yoga_wfs_child(const yoga_wfs_child& wfs_child);
  ~yoga_wfs_child();
};


class yoga_wfs {
 public:
  int                      device;
  int                      ndevices;
  int                      nhosts;
  int                      child_id;
  int                      slave_id;
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
  bool                     is_master;

  yoga_obj<cuFloatComplex> *d_camplipup;
  yoga_obj<cuFloatComplex> *d_camplifoc;
  yoga_obj<cuFloatComplex> *d_fttotim;
  yoga_obj<cuFloatComplex> *d_ftkernel;

  yoga_obj<float>          *d_pupil;
  yoga_obj<float>          *d_hrimg;
  yoga_obj<float>          *d_bincube;
  yoga_obj<float>          *d_binimg;
  yoga_obj<float>          *d_subsum;
  yoga_obj<float>          *d_offsets;
  yoga_obj<float>          *d_fluxPerSub;
  yoga_obj<float>          *d_sincar;
  yoga_obj<float>          *d_submask;
  yoga_obj<int>            *d_hrmap;

  yoga_obj<float>          *d_slopes;

  yoga_host_obj<float>     *image_telemetry;

  cufftHandle              pup_plan;
  cufftHandle              im_plan;

  // sh only
  yoga_obj<int>            *d_phasemap;
  yoga_obj<int>            *d_binmap;
  yoga_obj<int>            *d_validsubsx;  // nvalid
  yoga_obj<int>            *d_validsubsy;  // nvalid
  yoga_obj<int>            *d_istart;     // nxsub 
  yoga_obj<int>            *d_jstart;     // nxsub

  // pyramid only
  yoga_obj<float>          *d_psum;
  yoga_obj<cuFloatComplex> *d_phalfxy;
  yoga_obj<cuFloatComplex> *d_poffsets;

  yoga_host_obj<int>       *pyr_cx;
  yoga_host_obj<int>       *pyr_cy;

  yoga_source              *d_gs;

  yoga_streams 		   *streams;
  int  			   nstreams;

  yoga_context *current_context;

  vector<yoga_wfs *>       child_wfs;

 public:
  yoga_wfs(yoga_context *context, const char* type, long nxsub, long nvalid, long npix, long nphase, 
	   long nrebin, long nfft, long ntot, long npup,float pdiam,float nphotons, int lgs, int device);
  yoga_wfs(const yoga_wfs& wfs);
  ~yoga_wfs();

  int init_sh();
  int init_pyr();
  int wfs_initarrays(int *phasemap,int *hrmap, int *binmap,float *offsets, 
			     float *pupil, float *fluxPerSub, int *validsubsx, int *validsubsy, 
			     int *istart, int *jstart, cuFloatComplex *kernel);
  int wfs_initarrays(cuFloatComplex *halfxy,cuFloatComplex *offsets, float *focmask, 
		     float *pupil, int *cx, int *cy, float *sincar,int *phasemap,
		     int *validsubsx, int *validsubsy);
  int wfs_initgs(float xpos,float ypos,float lambda, float mag, long size,float noise, long seed);
  int load_kernels(float *lgskern);
  int sensor_trace(yoga_atmos *yatmos);
  int sensor_trace(yoga_dms *ydm, int rst);
  int comp_sh_generic();
  int comp_pyr_generic();
  int comp_sh_psf(cuFloatComplex *data_out, cuFloatComplex *data_inter, float *phase,
		  float *offsets, float *pupil, float scale,
		  int *istart,int *jstart,int *validx,int *validy,
		  int nphase, int screen_size, int nshim, int nssp,
		  cufftHandle plan, int device);
  int kernelconv_sh_psf(cuFloatComplex *data_out, cuFloatComplex *kernel, int nout, 
			cufftHandle plan, int device);
  int kernelconv_sh_psf(cuFloatComplex *data_out, cuFloatComplex *kernel, int nout, int nkernel,
			cufftHandle plan, int device);
  int comp_image_tele();
  int comp_image();
  int slopes_geom(int type, float *slopes);
  int slopes_geom(int type);
};

class yoga_sensors {
 public:
  int                 nsensors;
  vector<yoga_wfs *>  d_wfs;
     
 public:
  yoga_sensors(yoga_context *context, const char* type, int nwfs,long *nxsub,long *nvalid,long *npix,
	       long *nphase, long *nrebin,long *nfft, long *ntot, long npup, float *pdiam, float *nphot,
	       int *lgs, int device);
  ~yoga_sensors();

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
int fillbincube_async(yoga_streams *streams, float *bcube,cuFloatComplex  *hrimage, int *indxpix, int Nfft, 
		      int Npix, int Nrebin, int Nsub, int device);
int fillbinimg(float *bimage, float *bcube, int npix, int nsub, int Nsub, int *ivalid, int *jvalid, bool add, 
	       int device);
int fillbinimg_async(yoga_streams *streams, yoga_obj<float> *bimage, yoga_obj<float> *bcube, int npix, int nsub, 
		     int Nsub, int *ivalid, int *jvalid, bool add, int device);
int fillbinimg_async(yoga_host_obj<float> *image_telemetry, float *bimage, float *bcube, int npix,
		     int nsub, int Nsub, int *ivalid, int *jvalid,  int nim, bool add, int device);
int convolve(cuFloatComplex *d_odata,cuFloatComplex *d_idata,int N,int device);
int convolve_cube(cuFloatComplex *d_odata,cuFloatComplex *d_idata,int N,int n, int device);

// CUDA templates
// this is for cog
template <class T> void subap_reduce(int size, int threads, int blocks, T *d_idata, T *d_odata);
template <class T> void subap_reduce_async(int threads, int blocks, yoga_streams *streams, T *d_idata, T *d_odata);
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

void* comp_image_thread(void *thread_data);

#endif // _YOGA_WFS_H_
