// This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
//
// COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
// General Public License as published by the Free Software Foundation, either version 3 of the 
// License, or any later version.
//
// COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
// See the GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License along with COMPASS. 
// If not, see <https://www.gnu.org/licenses/>
//
//  Copyright (C) 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>

//! \file      sutra_controller.hpp
//! \ingroup   libsutra
//! \class     SutraController
//! \brief     this class provides the controller features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \version   5.5.0
//! \date      2022/01/24

#ifndef _SUTRA_CONTROLLER_H_
#define _SUTRA_CONTROLLER_H_

#include <carma_cublas.hpp>
#include <carma_host_obj.hpp>
#include <sutra_centroider.hpp>
#include <sutra_dm.hpp>
#include <sutra_utils.hpp>
#include <sutra_wfs.hpp>
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
  int32_t device;

  int32_t open_loop;
  Tcomp delay;
  Tcomp gain;
  float volt_min;
  float volt_max;
  int32_t nactus;
  int32_t nslopes;
  Tout val_max;
  Tcomp a;  // Coefficient for linear interpolation on command buffer to allow
            // non-integer delay
  Tcomp b;  // Coefficient for linear interpolation on command buffer to allow
            // non-integer delay
  Tcomp c;  // Coefficient for linear interpolation on command buffer to allow
            // non-integer delay
  vector<SutraDm *> d_dmseen;
  CarmaObj<Tcomp> *d_centroids;        // current centroids
  CarmaObj<Tcomp> *d_com;              // current command
  CarmaObj<Tcomp> *d_com_clipped;       // current command
  CarmaObj<Tout> *d_voltage;  // commands after perturbation and clipping
  CarmaObj<Tcomp> *d_com1;    // commands k-1
  vector<int32_t> centro_idx; // Centroider indices to handle
  std::deque<CarmaObj<Tcomp> *> d_circular_coms; //Circular buffer of commands for latency


  map<string, tuple<CarmaObj<Tcomp> *, int32_t, bool>> d_perturb_map;
  // perturbation command buffer

  cudaStream_t mainStream;

  // allocation of d_centroids and d_com
  SutraController(CarmaContext *context, int32_t nslope, int32_t nactu,
                   float delay, SutraDms *dms, int32_t *idx_dms, int32_t ndm,  int32_t *idx_centro, int32_t ncentro);
  virtual ~SutraController();

  virtual string get_type() = 0;

  //!!!! YOU MUST set d_centroids before calling it!!!!
  virtual int32_t comp_com() = 0;

  // It is better to have something like this (+protected d_centroids):
  // virtual int32_t comp_com (CarmaObj<T> *new_centroids)=0;
  // it would imply copy, but would be much safer

  inline int32_t nactu() { return this->nactus; }
  inline int32_t nslope() { return this->nslopes; }

  cublasHandle_t cublas_handle() { return current_context->get_cublas_handle(); }

  void init_voltage() {
    init_voltage_impl<Tcomp, Tout>(this->d_voltage, this->d_com_clipped);
  }

  int32_t set_centroids_ref(Tcomp *centroids_ref);
  int32_t add_perturb_voltage(string name, float *perturb, int32_t N);
  int32_t set_perturb_voltage(string name, float *perturb, int32_t N);
  int32_t remove_perturb_voltage(string name);
  int32_t reset_perturb_voltage();
  int32_t enable_perturb_voltage(string name);
  int32_t disable_perturb_voltage(string name);
  int32_t set_com(float *com, int32_t nElem);
  int32_t set_open_loop(int32_t open_loop_status, bool rst = true);
  int32_t clip_commands();
  int32_t comp_voltage();
  int32_t comp_latency();
  int32_t set_delay(float delay);
  int32_t set_volt_min(float volt_min);
  int32_t set_volt_max(float volt_max);
  int32_t set_val_max(float val_max);
  int32_t set_gain(float gain);
  int32_t reset_coms();

  // int32_t syevd_f(char meth, CarmaObj<T> *d_U,
  //             CarmaHostObj<T> *h_eingenvals);
  // int32_t invgen(CarmaObj<T> *d_mat, T cond, int32_t job);
  int32_t command_delay();
  int32_t add_perturb();

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
   * @return int32_t                             : error code
   */
  int32_t comp_polc(CarmaObj<Tcomp>& sk, CarmaObj<Tcomp>& iMat, CarmaObj<Tcomp>& ol_meas);

 protected:
  mutex comp_voltage_mutex;
};

template <typename Tin, typename Tout, std::enable_if_t<std::is_same<Tin, Tout>::value, bool> = true>
void convert_to_voltage(Tin *d_idata, Tout *d_odata, int32_t N, float volt_min, float volt_max,
                 uint16_t val_max, CarmaDevice *device, cudaStream_t stream){};

template <typename Tin, typename Tout, std::enable_if_t<!std::is_same<Tin, Tout>::value, bool> = true>
void convert_to_voltage(Tin *d_idata, Tout *d_odata, int32_t N, float volt_min, float volt_max,
                 uint16_t val_max, CarmaDevice *device, cudaStream_t stream);

int32_t shift_buf(float *d_data, int32_t offset, int32_t N, CarmaDevice *device);
int32_t fill_filtmat(float *filter, int32_t nactu, int32_t N, CarmaDevice *device);
int32_t TT_filt(float *mat, int32_t n, CarmaDevice *device);
int32_t fill_cmat(float *cmat, float *wtt, float *Mtt, int32_t nactu, int32_t nslopes,
              CarmaDevice *device);
int32_t do_statmat(float *statcov, int64_t dim, float *xpos, float *ypos, float norm,
               CarmaDevice *device);

template <class T>
int32_t get_pupphase(T *odata, float *idata, int32_t *indx_pup, int32_t Nphi,
                 CarmaDevice *device);

int32_t compute_Hcor_gpu(float *o_data, int32_t nrow, int32_t ncol, float Fs, float gmin,
                     float gmax, float delay, CarmaDevice *device);
int32_t absnormfft(cuFloatComplex *idata, float *odata, int32_t N, float norm,
               CarmaDevice *device);
int32_t adjust_csr_index(int32_t *rowind, int32_t *NNZ, int32_t *nact, int32_t nact_tot,
                     int32_t row_off, CarmaDevice *device);
template <typename T>
void pad_cmat(T *idata, int32_t m, int32_t n, T *odata, int32_t m2, int32_t n2,
              CarmaDevice *device);
#endif  // _SUTRA_CONTROLLER_H_
