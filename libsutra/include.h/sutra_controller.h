// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_controller.h
//! \ingroup   libsutra
//! \class     SutraController
//! \brief     this class provides the controller features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#ifndef _SUTRA_CONTROLLER_H_
#define _SUTRA_CONTROLLER_H_

#include <carma_cublas.h>
#include <carma_host_obj.h>
#include <sutra_centroider.h>
#include <sutra_dm.h>
#include <sutra_utils.h>
#include <sutra_wfs.h>
#include <mutex>
#include <tuple>
#include <deque>

using std::map;
using std::mutex;
using std::string;
using std::tuple;

template <typename Tcomp, typename Tout>
typename std::enable_if<std::is_same<Tcomp, Tout>::value, void>::type
init_voltage_impl(CarmaObj<Tout> *&volts, CarmaObj<Tcomp> *comClipped) {
  volts = comClipped;
};

template <typename Tcomp, typename Tout>
typename std::enable_if<!std::is_same<Tcomp, Tout>::value, void>::type
init_voltage_impl(CarmaObj<Tout> *&volts, CarmaObj<Tcomp> *comClipped) {
  volts = new CarmaObj<Tout>(comClipped->get_context(), comClipped->get_dims());
};

template <typename Tcomp, typename Tout>
class SutraController {
 public:
  CarmaContext *current_context;
  int device;

  int open_loop;
  Tcomp delay;
  Tcomp gain;
  float volt_min;
  float volt_max;
  int nactus;
  int nslopes;
  Tout val_max;
  Tcomp a;  // Coefficient for linear interpolation on command buffer to allow
            // non-integer delay
  Tcomp b;  // Coefficient for linear interpolation on command buffer to allow
            // non-integer delay
  Tcomp c;  // Coefficient for linear interpolation on command buffer to allow
            // non-integer delay
  vector<SutraDm *> d_dmseen;
  CarmaObj<Tcomp> *d_centroids;        // current centroids
  CarmaObj<Tcomp> *d_centroids_padded;  // current centroids
  CarmaObj<Tcomp> *d_com;              // current command
  CarmaObj<Tcomp> *d_com_padded;        // current command
  CarmaObj<Tcomp> *d_com_clipped;       // current command
  CarmaObj<Tout> *d_voltage;  // commands after perturbation and clipping
  CarmaObj<Tcomp> *d_com1;    // commands k-1
  vector<int> centro_idx; // Centroider indices to handle
  std::deque<CarmaObj<Tcomp> *> d_circular_coms; //Circular buffer of commands for latency


  map<string, tuple<CarmaObj<Tcomp> *, int, bool>> d_perturb_map;
  // perturbation command buffer

  cudaStream_t mainStream;

  // allocation of d_centroids and d_com
  SutraController(CarmaContext *context, int nslope, int nactu,
                   float delay, SutraDms *dms, int *idx_dms, int ndm,  int *idx_centro, int ncentro);
  virtual ~SutraController();

  virtual string get_type() = 0;

  //!!!! YOU MUST set d_centroids before calling it!!!!
  virtual int comp_com() = 0;

  // It is better to have something like this (+protected d_centroids):
  // virtual int comp_com (CarmaObj<T> *new_centroids)=0;
  // it would imply copy, but would be much safer

  inline int nactu() { return this->nactus; }
  inline int nslope() { return this->nslopes; }

  cublasHandle_t cublas_handle() { return current_context->get_cublas_handle(); }

  void init_voltage() {
    init_voltage_impl<Tcomp, Tout>(this->d_voltage, this->d_com_clipped);
  };

  int set_centroids_ref(Tcomp *centroids_ref);
  int add_perturb_voltage(string name, float *perturb, int N);
  int set_perturb_voltage(string name, float *perturb, int N);
  int remove_perturb_voltage(string name);
  int reset_perturb_voltage();
  int enable_perturb_voltage(string name);
  int disable_perturb_voltage(string name);
  int set_com(float *com, int nElem);
  int set_open_loop(int open_loop_status, bool rst = true);
  int clip_commands();
  int comp_voltage();
  int comp_latency();
  int set_delay(float delay);
  int set_volt_min(float volt_min);
  int set_volt_max(float volt_max);
  int set_val_max(float val_max);
  int set_gain(float gain);
  int reset_coms();

  // int syevd_f(char meth, CarmaObj<T> *d_U,
  //             CarmaHostObj<T> *h_eingenvals);
  // int invgen(CarmaObj<T> *d_mat, T cond, int job);
  int command_delay();
  int add_perturb();

  /**
   * @brief Compute the open loop measurements and effective commands
   *
   * eff_u = a * u_{k-1} + b * u_k
   * ol_meas = s_k - iMat * eff_u
   *
   * @param[in ] a       : Tcomp             :
   * @param[in ] uk_1    : CarmaObj<Tcomp>&  : commands at iteration k-1
   * @param[in ] b       : Tcomp             :
   * @param[in ] uk      : CarmaObj<Tcomp>&  : commands at iteration k
   * @param[in ] sk      : CarmaObj<Tcomp>&  : closed loop slopes at iteration k
   * @param[in ] iMat    : CarmaObj<Tcomp>&  : interaction matrix
   * @param[out] ol_meas : CarmaObj<Tcomp>&  : open loop measurements
   * @param[out] eff_u   : CarmaObj<Tcomp>&  : effective commands
   * @return int                             : error code
   */
  int comp_polc(CarmaObj<Tcomp>& sk, CarmaObj<Tcomp>& iMat, CarmaObj<Tcomp>& ol_meas);

 protected:
  mutex comp_voltage_mutex;
};

template <typename Tin, typename Tout, std::enable_if_t<std::is_same<Tin, Tout>::value, bool> = true>
void convert_to_voltage(Tin *d_idata, Tout *d_odata, int N, float volt_min, float volt_max,
                 uint16_t val_max, CarmaDevice *device, cudaStream_t stream){};

template <typename Tin, typename Tout, std::enable_if_t<!std::is_same<Tin, Tout>::value, bool> = true>
void convert_to_voltage(Tin *d_idata, Tout *d_odata, int N, float volt_min, float volt_max,
                 uint16_t val_max, CarmaDevice *device, cudaStream_t stream);

int shift_buf(float *d_data, int offset, int N, CarmaDevice *device);
int fill_filtmat(float *filter, int nactu, int N, CarmaDevice *device);
int TT_filt(float *mat, int n, CarmaDevice *device);
int fill_cmat(float *cmat, float *wtt, float *Mtt, int nactu, int nslopes,
              CarmaDevice *device);
int do_statmat(float *statcov, long dim, float *xpos, float *ypos, float norm,
               CarmaDevice *device);

template <class T>
int get_pupphase(T *odata, float *idata, int *indx_pup, int Nphi,
                 CarmaDevice *device);

int compute_Hcor_gpu(float *o_data, int nrow, int ncol, float Fs, float gmin,
                     float gmax, float delay, CarmaDevice *device);
int absnormfft(cuFloatComplex *idata, float *odata, int N, float norm,
               CarmaDevice *device);
int adjust_csr_index(int *rowind, int *NNZ, int *nact, int nact_tot,
                     int row_off, CarmaDevice *device);
#endif  // _SUTRA_CONTROLLER_H_
