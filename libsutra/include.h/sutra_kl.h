#ifndef _SUTRA_KL_H_
#define _SUTRA_KL_H_

#include <carma.h>
#include <carma_obj.h>
#include <carma_host_obj.h>

using namespace std;

class sutra_kl {
public:
	int device; // # device
	long dim;    // dim of final array
	long nr;     // # radial points
	long np;     // # of elements
	long nkl;    // # of functions in the basis

	carma_obj<float> *d_rabas; // the radial array of the basis
	carma_obj<float> *d_azbas; // the azimuthal array of the basis
	carma_host_obj<int> *h_ord;   // the radial orders of the basis
	carma_obj<int> *d_ord;   // the radial orders of the basis
	carma_obj<float> *d_cr;    //
	carma_obj<float> *d_cp;    //

  // Florian features
  carma_obj<float>          *d_covmat;
  carma_obj<float>          *d_filter;
  carma_obj<float>          *d_bas;
  carma_obj<float>          *d_evals;

  carma_context *current_context;   // the context in which it has been created 
public:
	sutra_kl(carma_context *context, long dim, long nr, long np, long nkl,
			int device);
	sutra_kl(const sutra_kl& kl);
	~sutra_kl();

	int do_compute(float alpha, float ampli, float *odata, int nkl, int size,
			int xoff, int yoff);
	int do_compute(float ampli, float *odata, int nkl, int size, int xoff,
			int yoff);
	int do_compute(float *odata, int nkl, int size, int xoff, int yoff);
	int do_combi(float *com, float *odata, int size, int xoff, int yoff);

// Florian features
  int get_flokl();
};
int getkl(float alpha, float ampli, float *d_odata, float *rabas, float *azbas,
		float *cr, float *cp, int nr, int np, int nx, int Nx, int xoff,
		int yoff);
int getkl(float ampli, float *d_odata, float *rabas, float *azbas, float *cr,
		float *cp, int nr, int np, int nx, int Nx, int xoff, int yoff);
int getkl(float *d_odata, float *rabas, float *azbas, float *cr, float *cp,
		int nr, int np, int nx, int Nx, int xoff, int yoff);
int combikl(float *com, int nkl, float *d_odata, float *rabas, int *h_ord,
		float *azbas, float *cr, float *cp, int nr, int np, int nx, int Nx,
		int xoff, int yoff);
int cget_flokl(long nkl, long dim, float *covmat, float *filter, float *bas);
//template <class T> void comp_kl(int threads, int blocks, T *d_idata, T *d_odata, int N);

#endif // _SUTRA_KL_H_
