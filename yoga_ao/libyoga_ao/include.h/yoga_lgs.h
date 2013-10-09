#ifndef _YOGA_LGS_H_
#define _YOGA_LGS_H_

#include <yoga.h>
#include <yoga_obj.h>
#include <yoga_ao_utils.h>

using namespace std;

class yoga_lgs {
 public:
  long                     nvalid;
  long                     npix;
  long                     nmaxhr;
  float                    hg;
  float                    h0;
  float                    deltah;
  float                    pixsize;
  long                     nprof;

  yoga_obj<float>          *d_doffaxis;
  yoga_obj<float>          *d_azimuth;
  yoga_obj<float>          *d_prof1d;
  yoga_obj<float>          *d_profcum;
  yoga_obj<cuFloatComplex> *d_prof2d;
  yoga_obj<float>          *d_beam;
  yoga_obj<cuFloatComplex> *d_ftbeam;
  yoga_obj<float>          *d_lgskern;
  yoga_obj<cuFloatComplex> *d_ftlgskern;
 
  yoga_context *current_context;
  /*
  cudaArray                *d_spotarray;
  cudaChannelFormatDesc    channelDesc;
  cudaMemcpy3DParms        copyParams;
  */

 public:
  yoga_lgs(yoga_context *context, long nvalid, long npix, long nmaxhr);
  yoga_lgs(const yoga_lgs& lgs);
  ~yoga_lgs();

  int lgs_init(int nprof, float hg, float h0, float deltah, float pixsie, float *doffaxis, float *prof1d,
	       float *profcum, float *beam, cuFloatComplex *ftbeam,float* azimuth);
  int load_prof(float *prof1d,float *profcum,float hg, float h0, float deltah);
  int lgs_update(int device);
  int lgs_makespot(int device, int nin);
  int load_kernels(float *lgskern, int device);

};


// General utilities
int interp_prof(cuFloatComplex *profout,float *prof1d,float *profcum, int npix, float *doffaxis, float hg, float pixsize, 
		 float h0, float deltah, int hmax, int Ntot, int device);
int times_ftbeam(cuFloatComplex *profout,cuFloatComplex *fbeam,int N, int Ntot, int device);
int rollbeamexp(float *imout, cuFloatComplex *iprof, float *beam,int N, int  Ntot,int device);
int lgs_rotate(cuFloatComplex *odata,float *idata, int width, int height, float *theta,float center,int Ntot,int device);
int rotate3d(cuFloatComplex *d_odata,cudaMemcpy3DParms copyParams, cudaArray *d_array, 
	     cudaChannelFormatDesc channelDesc, int width, int height, float *theta,float center, int Ntot,int device);

#endif // _YOGA_LGS_H_

