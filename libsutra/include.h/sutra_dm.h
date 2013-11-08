#ifndef _SUTRA_DM_H_
#define _SUTRA_DM_H_

#include <map>
#include <sutra_phase.h>
#include <sutra_kl.h>
#include <sutra_ao_utils.h>

typedef std::pair<std::string,float> type_screen;

using namespace std;

class sutra_dm {
 public:
  int                      device;
  string                   type;
  long                     ninflu;
  long                     influsize;
  long                     dim;
  float                    push4imat;

  sutra_phase               *d_shape;

  carma_obj<float>          *d_comm;

  carma_obj<float>          *d_influ; // if relevant

  // pzt 
  carma_obj<int>            *d_influpos;
  carma_obj<int>            *d_npoints;
  carma_obj<int>            *d_istart;
  carma_obj<int>            *d_xoff;
  carma_obj<int>            *d_yoff;

  //zernike
  carma_obj<float>          *d_coeffs;
  carma_obj<float>          *d_mask;
  carma_obj<float>          *d_zr;
  carma_obj<float>          *d_ztheta;

  sutra_kl      *d_kl;

  carma_context *current_context;

 public:
  sutra_dm(carma_context *context, const char* type, long dim, long ninflu,long influsize, long ninflupos, long n_npoints, float push4imat, int device);
  sutra_dm(const sutra_dm& dm);
  ~sutra_dm();

  int pzt_loadarrays(float *influ, int *influpos,int *npoints, int *istart, int *xoff, int *yoff);
  int kl_loadarrays(float *rabas, float *azbas,int *ord, float *cr, float *cp);
  int reset_shape();
  int comp_shape();
  int comp_shape(float *comm);
  int comp_oneactu(int nactu, float ampli);
};

class sutra_dms {
 public:
  int                         ndm;
  map<type_screen,sutra_dm *>  d_dms;    
     
 public:
  sutra_dms(int ndm);
  ~sutra_dms();

  int add_dm(carma_context *context, const char* type,float alt, long dim, long ninflu, long influsize, long ninflupos, long n_npoints, float push4imat, int device);
  int remove_dm(const char* type, float alt);
};

template <class T> void comp_dmshape(int threads, int blocks, T *d_idata, T *d_odata, int *pos, int *istart, int *npts, T *comm, unsigned int n, int N);
template <class T> void oneactu(int threads, int blocks, T *d_idata, T *d_odata, int nactu, T ampli, int *xoff, int *yoff, int dim_im, int dim_influ, int N);
template <class T> void oneactu(int threads, int blocks, T *d_idata, T *d_odata, int nactu, T ampli, int dim_im, int dim_influ, int N);
template <class T> void comp_fulldmshape(int threads, int blocks, T *d_idata, T *d_odata, int ninflu, int diminflu, T *comm, int N);
#endif // _SUTRA_DM_H_
