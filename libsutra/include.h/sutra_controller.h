#ifndef _sutra_controller_H_
#define _sutra_controller_H_

#include <carma_cublas.h>
#include <carma_host_obj.h>
#include <sutra_wfs.h>
#include <sutra_dm.h>
#include <sutra_centroider.h>
#include <sutra_ao_utils.h>

class sutra_controller {
public:

  //allocation of d_centroids and d_com
	sutra_controller(carma_context* context, int nslope, int nactu, float delay,
			sutra_dms *dms, char **type, float *alt, int ndm);
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
  set_centroids_ref(float *centroids_ref);
	int
  get_centroids_ref(float *centroids_ref);
	int
  set_perturbcom(float *perturb, int N);
  int
  set_openloop(int open_loop_status);
  void
  clip_voltage(float min, float max);
  int
  comp_voltage();
  int
  syevd_f(char meth, carma_obj<float> *d_U, carma_host_obj<float> *h_eingenvals);
  int
  invgen(carma_obj<float> *d_mat, float cond, int job);
  int
  command_delay();

public:
//I would propose to make them protected (+ proper
//set of fuctions). It could make life easier!
//But we should discuss it
  int cpt_pertu;
  int open_loop;
  float delay;
  float a; // Coefficient for linear interpolation on command buffer to allow non-integer delay
  float b; // Coefficient for linear interpolation on command buffer to allow non-integer delay
  float c; // Coefficient for linear interpolation on command buffer to allow non-integer delay
  vector<sutra_dm *> d_dmseen;
  carma_obj<float> *d_subsum; // current flux
	carma_obj<float> *d_centroids; // current centroids
	carma_obj<float> *d_centroids_ref; // ref centroids
  carma_obj<float> *d_com; // current command
  carma_obj<float> *d_perturb; // perturbation command buffer
  carma_obj<float> *d_voltage; // commands sent to mirror
  carma_obj<float> *d_com1; // commands k-1
  carma_obj<float> *d_com2; // commands k-2

  carma_streams *streams;

protected:
  int device;
  carma_context *current_context;

};

int
shift_buf(float *d_data, int offset, int N, carma_device *device);
int
mult_vect(float *d_data, float *scale, int N, carma_device *device);
int
mult_vect(float *d_data, float *scale, float gain, int N, carma_device *device);
int
mult_vect(float *d_data, float gain, int N, carma_device *device);
int
mult_int(float *o_data, float *i_data, float *scale, float gain, int N,
    carma_device *device);
int
mult_int(float *o_data, float *i_data, float *scale, float gain, int N,
    carma_device *device, carma_streams *streams);
int
mult_int(float *o_data, float *i_data, float gain, int N,carma_device *device);
int
fill_filtmat(float *filter, int nactu, int N, carma_device *device);
int
TT_filt(float *mat, int n, carma_device *device);
int
fill_cmat(float *cmat, float *wtt, float *Mtt, int nactu, int nslopes, carma_device *device);
int
do_statmat(float *statcov,long dim, float *xpos, float *ypos, float norm, carma_device *device);
int
add_md(float *o_matrix, float *i_matrix, float *i_vector, int N,carma_device *device);

template<class T>
int
get_pupphase(T *odata, float *idata, int *indx_pup, int Nphi, carma_device *device);

int
compute_Hcor_gpu(float *o_data, int nrow, int ncol, float Fs, float gmin, float gmax, float delay, carma_device *device);
int
absnormfft(cuFloatComplex *idata, float *odata, int N, float norm, carma_device *device);
int
adjust_csr_index(int *rowind, int *NNZ, int *nact, int nact_tot, int row_off, carma_device *device);
#endif // _sutra_controller_H_
