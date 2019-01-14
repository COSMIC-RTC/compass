#ifndef _sutra_controller_H_
#define _sutra_controller_H_

#include <carma_cublas.h>
#include <carma_host_obj.h>
#include <sutra_centroider.h>
#include <sutra_dm.h>
#include <sutra_utils.h>
#include <sutra_wfs.h>
#include <mutex>
#include <tuple>

using std::map;
using std::mutex;
using std::string;
using std::tuple;

class sutra_controller {
 public:
  carma_context *current_context;
  int device;

  // allocation of d_centroids and d_com
  sutra_controller(carma_context *context, int nvalid, int nslope, int nactu,
                   float delay, sutra_dms *dms, int *idx_dms, int ndm);
  virtual ~sutra_controller();

  virtual string get_type() = 0;

  //!!!! YOU MUST set d_centroids before calling it!!!!
  virtual int comp_com() = 0;

  // It is better to have something like this (+protected d_centroids):
  // virtual int comp_com (carma_obj<float> *new_centroids)=0;
  // it would imply copy, but would be much safer

  inline int nactu() { return d_com->getDims(1); }
  inline int nslope() { return d_centroids->getDims(1); }

  cublasHandle_t cublas_handle() { return current_context->get_cublasHandle(); }

  int set_centroids_ref(float *centroids_ref);
  int get_centroids_ref(float *centroids_ref);
  int add_perturb_voltage(string name, float *perturb, int N);
  int remove_perturb_voltage(string name);
  int reset_perturb_voltage();
  int enable_perturb_voltage(string name);
  int disable_perturb_voltage(string name);
  int set_com(float *com, int nElem);
  int set_openloop(int open_loop_status, bool rst = true);
  void clip_voltage(float min, float max);
  int comp_voltage();
  int syevd_f(char meth, carma_obj<float> *d_U,
              carma_host_obj<float> *h_eingenvals);
  int invgen(carma_obj<float> *d_mat, float cond, int job);
  int command_delay();
  int add_perturb();

 public:
  // I would propose to make them protected (+ proper
  // set of fuctions). It could make life easier!
  // But we should discuss it
  int open_loop;
  float delay;
  float a;  // Coefficient for linear interpolation on command buffer to allow
            // non-integer delay
  float b;  // Coefficient for linear interpolation on command buffer to allow
            // non-integer delay
  float c;  // Coefficient for linear interpolation on command buffer to allow
            // non-integer delay
  vector<sutra_dm *> d_dmseen;
  carma_obj<float> *d_subsum;         // current flux
  carma_obj<float> *d_centroids;      // current centroids
  carma_obj<float> *d_centroids_ref;  // ref centroids
  carma_obj<float> *d_com;            // current command
  carma_obj<float> *d_voltage;        // commands sent to mirror
  carma_obj<float> *d_com1;           // commands k-1
  carma_obj<float> *d_com2;           // commands k-2

  map<string, tuple<carma_obj<float> *, int, bool>> d_perturb_map;
  // perturbation command buffer

  carma_streams *streams;

 protected:
  mutex comp_voltage_mutex;
};

int shift_buf(float *d_data, int offset, int N, carma_device *device);
int fill_filtmat(float *filter, int nactu, int N, carma_device *device);
int TT_filt(float *mat, int n, carma_device *device);
int fill_cmat(float *cmat, float *wtt, float *Mtt, int nactu, int nslopes,
              carma_device *device);
int do_statmat(float *statcov, long dim, float *xpos, float *ypos, float norm,
               carma_device *device);

template <class T>
int get_pupphase(T *odata, float *idata, int *indx_pup, int Nphi,
                 carma_device *device);

int compute_Hcor_gpu(float *o_data, int nrow, int ncol, float Fs, float gmin,
                     float gmax, float delay, carma_device *device);
int absnormfft(cuFloatComplex *idata, float *odata, int N, float norm,
               carma_device *device);
int adjust_csr_index(int *rowind, int *NNZ, int *nact, int nact_tot,
                     int row_off, carma_device *device);
#endif  // _sutra_controller_H_
