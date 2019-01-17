#ifndef _SUTRA_WFS_H_
#define _SUTRA_WFS_H_

#include <carma_utils.h>
#include <sutra_lgs.h>
#include <sutra_phase.h>
#include <sutra_target.h>
#include <sutra_telemetry.h>
#include <sutra_telescope.h>
#include <sutra_utils.h>

#include <map>
#include <vector>
//#include <sutra_slopes.h>

using std::map;
using std::string;

class sutra_wfs {
 public:
  int device;
  string type;
  long nxsub;
  long nvalid;
  long npix;
  long nrebin;
  long nfft;
  long ntot;
  long npup;
  long nphase;
  long nmaxhr;
  long nffthr;
  float subapd;
  float nphot;
  float nphot4imat;
  float noise;
  bool lgs;
  bool kernconv;
  bool roket;
  bool is_low_order;

  cufftHandle *campli_plan;
  cufftHandle *fttotim_plan;
  carma_obj<cuFloatComplex> *d_ftkernel;
  carma_obj<cuFloatComplex> *d_camplipup;
  carma_obj<cuFloatComplex> *d_camplifoc;
  carma_obj<cuFloatComplex> *d_fttotim;

  carma_obj<float> *d_pupil;
  carma_obj<float> *d_bincube;
  carma_obj<float> *d_binimg;
  carma_obj<float> *d_binimg_notnoisy;
  carma_obj<float> *d_intensities;
  carma_obj<float> *d_offsets;
  carma_obj<float> *d_fluxPerSub;
  carma_obj<float> *d_sincar;
  carma_obj<int> *d_hrmap;

  carma_obj<float> *d_slopes;

  carma_host_obj<float> *image_telemetry;

  sutra_source *d_gs;
  std::vector<carma_obj<float> *> d_pupil_ngpu;

  carma_streams *streams;
  int nstreams;

  carma_obj<int> *d_phasemap;
  carma_obj<int> *d_validsubsx;  // nvalid
  carma_obj<int> *d_validsubsy;  // nvalid

  carma_context *current_context;

  /// MPI stuff
  int offset;
  int nvalid_tot;
  int rank;
  int *displ_bincube;
  int *count_bincube;

 public:
  virtual ~sutra_wfs(){};

  int wfs_initgs(carma_obj<float> *d_lgskern,
                 carma_obj<cuFloatComplex> *d_ftlgskern,
                 map<vector<int>, cufftHandle *> ftlgskern_plans, float xpos,
                 float ypos, float lambda, float mag, float zerop, long size,
                 float noise, long seed, float G, float thetaML, float dx,
                 float dy);
  int set_pupil(float *pupil);
  int set_binimg(float *binimg, int nElem);
  int load_kernels(float *lgskern);
  int sensor_trace(sutra_atmos *yatmos);
  int sensor_trace(sutra_dms *ydm, int rst);
  int sensor_trace(sutra_atmos *atmos, sutra_dms *ydms);
  int sensor_trace(int rst);
  int slopes_geom(float *slopes, int type = 0);
  int slopes_geom(int type = 0);

  virtual int fill_binimage(int async) = 0;
  virtual int comp_image(bool noise = true) = 0;

  virtual int define_mpi_rank(int rank, int size) = 0;
  virtual int allocate_buffers(
      map<vector<int>, cufftHandle *> campli_plans,
      map<vector<int>, cufftHandle *> fttotim_plans) = 0;
  int set_noise(float noise, long seed);

 protected:
  virtual int comp_generic() = 0;
  sutra_wfs(carma_context *context, sutra_telescope *d_tel,
            carma_obj<cuFloatComplex> *d_camplipup,
            carma_obj<cuFloatComplex> *d_camplifoc,
            carma_obj<cuFloatComplex> *d_fttotim, string type, long nxsub,
            long nvalid, long npix, long nphase, long nrebin, long nfft,
            long ntot, long npup, float pdiam, float nphotons, float nphot4imat,
            int lgs, bool is_low_order, bool roket, int device);
};

// General utilities
int fillcamplipup(cuFloatComplex *amplipup, float *phase, float *offset,
                  float *mask, float scale, int *istart, int *jstart,
                  int *ivalid, int *jvalid, int nphase, int npup, int Nfft,
                  int Ntot, carma_device *device, int offset_phase);
int indexfill(cuFloatComplex *d_odata, cuFloatComplex *d_idata, int *indx,
              int nx, int Nx, int N, carma_device *device);
int fillbincube(float *bcube, cuFloatComplex *hrimage, int *indxpix, int Nfft,
                int Npix, int Nrebin, int Nsub, carma_device *device);
int fillbincube_async(carma_streams *streams, float *bcube,
                      cuFloatComplex *hrimage, int *indxpix, int Nfft, int Npix,
                      int Nrebin, int Nsub, carma_device *device);
int fillbinimg(float *bimage, float *bcube, int npix, int nsub, int Nsub,
               int *ivalid, int *jvalid, bool add, carma_device *device);
int fillbinimg_async(carma_streams *streams, carma_obj<float> *bimage,
                     carma_obj<float> *bcube, int npix, int nsub, int Nsub,
                     int *ivalid, int *jvalid, bool add, carma_device *device);
int fillbinimg_async(carma_host_obj<float> *image_telemetry, float *bimage,
                     float *bcube, int npix, int nsub, int Nsub, int *ivalid,
                     int *jvalid, int nim, bool add, carma_device *device);
int convolve_cube(cuFloatComplex *d_odata, cuFloatComplex *d_idata, int N,
                  int n, carma_device *device);

// CUDA templates
// this is for cog
template <class T>
void subap_reduce(int size, int threads, int blocks, T *d_idata, T *d_odata,
                  carma_device *device);

template <class T>
void subap_reduce_async(int threads, int blocks, carma_streams *streams,
                        T *d_idata, T *d_odata);
// this is for tcog
template <class T>
void subap_reduce(int size, int threads, int blocks, T *d_idata, T *d_odata,
                  T thresh, carma_device *device);
// this is for tcog_new
template <class T>
void subap_reduce_new(int size, int threads, int blocks, T *d_idata, T *d_odata,
                      T thresh, carma_device *device);
// this is for wcog
template <class T>
void subap_reduce(int size, int threads, int blocks, T *d_idata, T *d_odata,
                  T *weights, carma_device *device);

template <class T>
void phase_reduce(int threads, int blocks, T *d_idata, T *d_odata, int *indx,
                  T alpha);

template <class T>
void phase_derive(int size, int threads, int blocks, int n, T *d_idata,
                  T *d_odata, int *indx, T *mask, T alpha, float *fluxPerSub);

template <class Tout, class Tin>
void pyr_getpup(Tout *d_odata, Tin *d_idata, Tout *d_offsets, Tin *d_pup,
                int np, float lambda, carma_device *device);

template <class Tout, class Tin>
void pyr_getpup(Tout *d_odata, Tin *d_idata, Tin *d_pup, int np, int N,
                float lambda, float cx, float cy, carma_device *device);

template <class T>
void pyr_rollmod(T *d_odata, T *d_idata, T *d_mask, float cx, float cy, int np,
                 int ns, carma_device *device);

template <class T>
void pyr_fillbin(T *d_odata, T *d_idata, int nrebin, int np, int ns, int nim,
                 carma_device *device);

template <class T>
int pyr_fillbinimg(T *bimage, const T *bcube, const int nxsub, const bool add,
                   carma_device *device);

template <class T>
int pyr_fillbinimg(T *oimage, const T *image, const int n, const int N,
                   const int rebin, const bool add, carma_device *device);

template <class Tin, class Tout>
void pyr_abs2(Tout *d_odata, Tin *d_idata, Tout fact, int ns, int nim,
              carma_device *device);

template <class Tout, class Tin>
void pyr_submask(Tout *d_odata, Tin *d_mask, int n, carma_device *device);

template <class T>
void pyr_submaskpyr(T *d_odata, T *d_mask, int n, carma_device *device);

template <class T>
void pyr_submaskpyr(T *d_odata, float *d_mask, int n, carma_device *device);

template <class Tout, class Tin>
void pyr_abs(Tout *d_odata, Tin *d_idata, int ns, int nim,
             carma_device *device);

template <class Tout, class Tin>
void pyr_submask3d(Tout *d_odata, Tin *d_mask, int n, int nim,
                   carma_device *device);

template <class T>
void pyr_intensities(T *d_odata, T *d_idata, int *subindx, int *subindy, int ns,
                     int nvalid, int nim, carma_device *device);

template <class T>
void pyr_intensities(T *d_odata, T *d_idata, int *subindx, int *subindy, int ns,
                     int nvalid, carma_device *device);

template <class T>
void pyr_fact(T *d_data, T fact, int n, int nim, carma_device *device);

void pyr_fact(cuFloatComplex *d_data, float fact, int n, int nim,
              carma_device *device);
void pyr_fact(float *d_data, float fact1, float *fact2, int n, int nim,
              carma_device *device);

template <class Tin, class Tout>
void roof_abs2(Tout *d_odata, Tin *d_idata, Tout fact, int ns, int nim,
               carma_device *device);

template <class T>
void roof_intensities(T *d_odata, T *d_idata, int *subindx, int *subindy,
                      int ns, int nvalid, int nim, carma_device *device);

template <class T>
void roof_rollmod(T *d_odata, T *d_idata, T *d_mask, float cx, float cy, int np,
                  int ns, carma_device *device);

template <class T>
void roof_fillbin(T *d_odata, T *d_idata, int nrebin, int np, int ns, int nim,
                  carma_device *device);

void copyImginBinimg(float *binimg, int *validsubsx, int *validsubsy, int Nb,
                     float *img, int *validx, int *validy, int Nim, int Npix,
                     carma_device *device);

#endif  // _SUTRA_WFS_H_
