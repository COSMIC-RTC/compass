#ifndef _sutra_controller_H_
#define _sutra_controller_H_

#include <carma_cublas.h>
#include <carma_host_obj.h>
#include <sutra_wfs.h>
#include <sutra_dm.h>
#include <sutra_centroider.h>
#include <sutra_ao_utils.h>

using namespace std;

class sutra_controller {
public:

  //allocation of d_centroids and d_com
  sutra_controller(carma_context* context, int nslope, int nactu);
  virtual
  ~sutra_controller();

  virtual string
  get_type()=0;

  //!!!! YOU MUST set d_centroids before calling it!!!!
  virtual int
  comp_com()=0;

  //It is better to have something like this (+protected d_centroids):
  //virtual int comp_com (carma_obj<float> *new_centroids)=0;
  //it would imply copy, but would be much safer

  inline int nactu() {
    return d_com->getDims(1);
  }
  inline int nslope() {
    return d_centroids->getDims(1);
  }

  cublasHandle_t cublas_handle() {
    return current_context->get_cublasHandle();
  }

  int
  syevd_f(char meth, carma_obj<float> *d_U, carma_host_obj<float> *h_eingenvals);
  int
  invgen(carma_obj<float> *d_mat, float cond, int job);

public:
//I would propose to make them protected (+ proper
//set of fuctions). It could make life easier!
//But we should discuss it
  carma_obj<float> *d_centroids; // current centroids
  carma_obj<float> *d_com; // current command

  carma_streams *streams;

protected:
  int device;
  carma_context *current_context;

};

int
shift_buf(float *d_data, int offset, int N, int device);
int
mult_vect(float *d_data, float *scale, int N, int device);
int
mult_vect(float *d_data, float gain, int N, int device);
int
mult_int(float *o_data, float *i_data, float *scale, float gain, int N,
    int device);
int
mult_int(float *o_data, float *i_data, float *scale, float gain, int N,
    int device, carma_streams *streams);
int
mult_int(float *o_data, float *i_data, float gain, int N,int device);
int
fill_filtmat(float *filter, int nactu, int N, int device);
int
TT_filt(float *mat, int n, int device);
int
fill_cmat(float *cmat, float *wtt, float *Mtt, int nactu, int nslopes, int device);
int
do_statmat(float *statcov,long dim, float *xpos, float *ypos, float norm, int device);
int
add_md(float *o_matrix, float *i_matrix, float *i_vector, int N,int device);
int
floattodouble(float *idata, double *odata, int N, int device);
int
doubletofloat(double *idata, float *odata, int N, int device);
int
get_pupphase(float *odata, float *idata, int *indx_pup, int Nphi, int device);
int
compute_Hcor_gpu(float *o_data, int nrow, int ncol, float Fs, float gmin, float gmax, int delay, int device);
int
absnormfft(cuFloatComplex *idata, float *odata, int N, float norm, int device);
#endif // _sutra_controller_H_
