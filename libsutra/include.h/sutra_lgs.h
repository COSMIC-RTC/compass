#ifndef _SUTRA_LGS_H_
#define _SUTRA_LGS_H_

#include <carma.h>
#include <carma_obj.h>
#include <sutra_ao_utils.h>

using namespace std;

class sutra_sensors;
class sutra_lgs {
public:
  long nvalid;
  long npix;
  long nmaxhr;
  float hg;
  float h0;
  float deltah;
  float pixsize;
  long nprof;

  cufftHandle *ftlgskern_plan;
  carma_obj<float> *d_doffaxis;
  carma_obj<float> *d_azimuth;
  carma_obj<float> *d_prof1d;
  carma_obj<float> *d_profcum;
  carma_obj<cuFloatComplex> *d_prof2d;
  carma_obj<float> *d_beam;
  carma_obj<cuFloatComplex> *d_ftbeam;
  carma_obj<float> *d_lgskern;
  carma_obj<cuFloatComplex> *d_ftlgskern;

  carma_context *current_context;
  /*
   cudaArray                *d_spotarray;
   cudaChannelFormatDesc    channelDesc;
   cudaMemcpy3DParms        copyParams;
   */

public:
  sutra_lgs(carma_context *context, sutra_sensors *sensors, long nvalid, long npix, long nmaxhr);
  sutra_lgs(const sutra_lgs& lgs);
  ~sutra_lgs();

  int
  lgs_init(int nprof, float hg, float h0, float deltah, float pixsie,
      float *doffaxis, float *prof1d, float *profcum, float *beam,
      cuFloatComplex *ftbeam, float* azimuth);
  int
  load_prof(float *prof1d, float *profcum, float hg, float h0, float deltah);
  int
  lgs_update(carma_device *device);
  int
  lgs_makespot(carma_device *device, int nin);
  int
  load_kernels(float *lgskern, carma_device *device);

};

// General utilities
int
interp_prof(cuFloatComplex *profout, float *prof1d, float *profcum, int npix,
    float *doffaxis, float hg, float pixsize, float h0, float deltah, int hmax,
    int Ntot, carma_device *device);
int
times_ftbeam(cuFloatComplex *profout, cuFloatComplex *fbeam, int N, int Ntot,
    carma_device *device);
int
rollbeamexp(float *imout, cuFloatComplex *iprof, float *beam, int N, int Ntot,
    carma_device *device);
int
lgs_rotate(cuFloatComplex *odata, float *idata, int width, int height,
    float *theta, float center, int Ntot, carma_device *device);
int
rotate3d(cuFloatComplex *d_odata, cudaMemcpy3DParms copyParams,
    cudaArray *d_array, cudaChannelFormatDesc channelDesc, int width,
    int height, float *theta, float center, int Ntot, carma_device *device);

#endif // _SUTRA_LGS_H_
