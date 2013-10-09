#ifndef _YOGA_KL_H_
#define _YOGA_KL_H_

#include <yoga.h>
#include <yoga_obj.h>
#include <yoga_host_obj.h>

using namespace std;

class yoga_kl {
 public:
  int                      device; // # device
  long                     dim;    // dim of final array
  long                     nr;     // # radial points
  long                     np;     // # of elements
  long                     nkl;    // # of functions in the basis

  yoga_obj<float>          *d_rabas; // the radial array of the basis
  yoga_obj<float>          *d_azbas; // the azimuthal array of the basis
  yoga_host_obj<int>       *h_ord;   // the radial orders of the basis
  yoga_obj<int>            *d_ord;   // the radial orders of the basis
  yoga_obj<float>          *d_cr;    //
  yoga_obj<float>          *d_cp;    //


  yoga_context *current_context;   // the context in which it has been created 

 public:
  yoga_kl(yoga_context *context, long dim, long nr, long np, long nkl, int device);
  yoga_kl(const yoga_kl& kl);
  ~yoga_kl();

  int do_compute(float alpha, float ampli, float *odata, int nkl, int size, int xoff, int yoff);
  int do_compute(float ampli, float *odata, int nkl, int size, int xoff, int yoff);
  int do_compute(float *odata, int nkl, int size, int xoff, int yoff);
  int do_combi(float *com, float *odata, int size, int xoff, int yoff);
};
int getkl(float alpha, float ampli, float *d_odata,float *rabas, float *azbas, float *cr, float *cp, int nr, int np, int nx, int Nx, int xoff, int yoff);
int getkl(float ampli, float *d_odata,float *rabas, float *azbas, float *cr, float *cp, int nr, int np, int nx, int Nx, int xoff, int yoff);
int getkl(float *d_odata,float *rabas, float *azbas, float *cr, float *cp, int nr, int np, int nx, int Nx, int xoff, int yoff);
int combikl(float *com, int nkl, float *d_odata,float *rabas, int *h_ord, float *azbas, float *cr, float *cp, int nr, int np, int nx, int Nx, int xoff, int yoff);
//template <class T> void comp_kl(int threads, int blocks, T *d_idata, T *d_odata, int N);

#endif // _YOGA_KL_H_
