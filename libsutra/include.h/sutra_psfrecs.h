#ifndef _SUTRA_PSFRECS_H_
#define _SUTRA_PSFRECS_H_

#include <carma.h>
#include <carma_obj.h>
#include <carma_sparse_obj.h>
#include <carma_host_obj.h>

class sutra_psfrecs {
  public:
    carma_context *current_context;
    int device;

    int nactus; // number of actuators
    int niter; // number of iterations
    int nmodes; // number of modes
    //PSF reconstruction from roket data
    carma_obj<float> *d_err;
    carma_obj<cuFloatComplex> *d_amplipup;
    carma_obj<float> *d_psf;
    carma_obj<float> *d_phase;
    carma_obj<int> *d_wherephase;
    carma_sparse_obj<float> *d_IF;
    carma_obj<float> *d_TT;
    float scale;
    int size;
    int Npts;
    // PSF reconstruction with Vii functions
    carma_obj<float> *d_term1;
    carma_obj<float> *d_term2;
    carma_obj<float> *d_otftel;
    carma_obj<float> *d_otfVii;
    carma_obj<float> *d_mask;
    carma_obj<float> *d_eigenvals;
    carma_obj<float> *d_Btt;
    carma_obj<float> *d_covmodes;
    carma_host_obj<float> *h_eigenvals;
    carma_obj<cuFloatComplex> *d_newmodek;
    carma_obj<cuFloatComplex> *d_Dphi;
    carma_obj<cuFloatComplex> *d_pupfft;


  public:
    sutra_psfrecs(carma_context *context, int device, char *type, int nactus,
                int nmodes, int niter, float *IFvalue, int *IFrowind, int *IFcolind, int IFnz,
                float *d_TT, float *pupil, int size, int Npts, float scale, float *Btt, float *covmodes);
    ~sutra_psfrecs();

    int
    psf_rec_roket(float *err);
    int
    psf_rec_Vii();
};

int
fill_amplipup(cuFloatComplex *amplipup, float *phase, int *wherephase,
                float scale, int Npts, int nx, int Nx, int puponly, carma_device *device);
int
cumulpsf(float *d_odata, cuFloatComplex *d_idata, int N, carma_device *device);
int
abs2complex(cuFloatComplex *d_odata, cuFloatComplex *d_idata, int N, carma_device *device);
int
real(float *d_odata, cuFloatComplex *d_idata, int N, carma_device *device);
int
fill_mask(float *d_odata, float *d_idata, int N, int norm, carma_device *device);
int
modulus2(float *d_odata, cuFloatComplex *d_idata, int N, carma_device *device);
int
pow2(cuFloatComplex *d_odata, cuFloatComplex *d_idata, int N, carma_device *device);
int
fill_term1(float *d_odata, cuFloatComplex *d_idata, cuFloatComplex *d_pupfft, int N, carma_device *device);
int
add2Dphi(cuFloatComplex *d_odata, float *d_term1, float *d_term2, float e, int N, carma_device *device);
int
computeOTFvii(float *d_otfVii, cuFloatComplex * d_Dphi, float *d_otftel, float *d_mask,
                    float scale, int N, carma_device *device);
int
ifftscale(cuFloatComplex *d_odata, float scale, int N, carma_device *device);

#endif
