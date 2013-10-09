#ifndef _YOGA_DM_H_
#define _YOGA_DM_H_

#include <map>
#include <yoga_phase.h>
#include <yoga_kl.h>
#include <yoga_ao_utils.h>

typedef std::pair<std::string,float> type_screen;

using namespace std;

class yoga_dm {
 public:
  int                      device;
  string                   type;
  long                     ninflu;
  long                     influsize;
  long                     dim;
  float                    push4imat;

  yoga_phase               *d_shape;

  yoga_obj<float>          *d_comm;

  yoga_obj<float>          *d_influ; // if relevant

  // pzt 
  yoga_obj<int>            *d_influpos;
  yoga_obj<int>            *d_npoints;
  yoga_obj<int>            *d_istart;
  yoga_obj<int>            *d_xoff;
  yoga_obj<int>            *d_yoff;

  //zernike
  yoga_obj<float>          *d_coeffs;
  yoga_obj<float>          *d_mask;
  yoga_obj<float>          *d_zr;
  yoga_obj<float>          *d_ztheta;

  yoga_kl      *d_kl;

  yoga_context *current_context;

 public:
  yoga_dm(yoga_context *context, const char* type, long dim, long ninflu,long influsize, long ninflupos, long n_npoints, float push4imat, int device);
  yoga_dm(const yoga_dm& dm);
  ~yoga_dm();

  int pzt_loadarrays(float *influ, int *influpos,int *npoints, int *istart, int *xoff, int *yoff);
  int kl_loadarrays(float *rabas, float *azbas,int *ord, float *cr, float *cp);
  int reset_shape();
  int comp_shape();
  int comp_shape(float *comm);
  int comp_oneactu(int nactu, float ampli);
};

class yoga_dms {
 public:
  int                         ndm;
  map<type_screen,yoga_dm *>  d_dms;    
     
 public:
  yoga_dms(int ndm);
  ~yoga_dms();

  int add_dm(yoga_context *context, const char* type,float alt, long dim, long ninflu, long influsize, long ninflupos, long n_npoints, float push4imat, int device);
  int remove_dm(const char* type, float alt);
};

template <class T> void comp_dmshape(int threads, int blocks, T *d_idata, T *d_odata, int *pos, int *istart, int *npts, T *comm, unsigned int n, int N);
template <class T> void oneactu(int threads, int blocks, T *d_idata, T *d_odata, int nactu, T ampli, int *xoff, int *yoff, int dim_im, int dim_influ, int N);
template <class T> void oneactu(int threads, int blocks, T *d_idata, T *d_odata, int nactu, T ampli, int dim_im, int dim_influ, int N);
template <class T> void comp_fulldmshape(int threads, int blocks, T *d_idata, T *d_odata, int ninflu, int diminflu, T *comm, int N);
#endif // _YOGA_DM_H_
