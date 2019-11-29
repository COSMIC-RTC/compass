// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2019 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.
//  Distributed under GNU - LGPL
//
//  COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser 
//  General Public License as published by the Free Software Foundation, either version 3 of the License, 
//  or any later version.
//
//  COMPASS: End-to-end AO simulation tool using GPU acceleration 
//  The COMPASS platform was designed to meet the need of high-performance for the simulation of AO systems. 
//  
//  The final product includes a software package for simulating all the critical subcomponents of AO, 
//  particularly in the context of the ELT and a real-time core based on several control approaches, 
//  with performances consistent with its integration into an instrument. Taking advantage of the specific 
//  hardware architecture of the GPU, the COMPASS tool allows to achieve adequate execution speeds to
//  conduct large simulation campaigns called to the ELT. 
//  
//  The COMPASS platform can be used to carry a wide variety of simulations to both testspecific components 
//  of AO of the E-ELT (such as wavefront analysis device with a pyramid or elongated Laser star), and 
//  various systems configurations such as multi-conjugate AO.
//
//  COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the 
//  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
//  See the GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License along with COMPASS. 
//  If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>.
// -----------------------------------------------------------------------------

//! \file      sutra_controller.h
//! \ingroup   libsutra
//! \class     sutra_controller
//! \brief     this class provides the controller features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.3.2
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

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
init_voltage_impl(carma_obj<Tout> *&volts, carma_obj<Tcomp> *comClipped) {
  volts = comClipped;
};

template <typename Tcomp, typename Tout>
typename std::enable_if<!std::is_same<Tcomp, Tout>::value, void>::type
init_voltage_impl(carma_obj<Tout> *&volts, carma_obj<Tcomp> *comClipped) {
  volts = new carma_obj<Tout>(comClipped->getContext(), comClipped->getDims());
};

template <typename Tcomp, typename Tout>
class sutra_controller {
 public:
  carma_context *current_context;
  int device;

  int open_loop;
  Tcomp delay;
  Tcomp gain;
  float Vmin;
  float Vmax;
  int nactus;
  int nslopes;
  Tout valMax;
  Tcomp a;  // Coefficient for linear interpolation on command buffer to allow
            // non-integer delay
  Tcomp b;  // Coefficient for linear interpolation on command buffer to allow
            // non-integer delay
  Tcomp c;  // Coefficient for linear interpolation on command buffer to allow
            // non-integer delay
  vector<sutra_dm *> d_dmseen;
  carma_obj<Tcomp> *d_centroids;        // current centroids
  carma_obj<Tcomp> *d_centroidsPadded;  // current centroids
  carma_obj<Tcomp> *d_com;              // current command
  carma_obj<Tcomp> *d_comPadded;        // current command
  carma_obj<Tcomp> *d_comClipped;       // current command
  carma_obj<Tout> *d_voltage;  // commands after perturbation and clipping
  carma_obj<Tcomp> *d_com1;    // commands k-1
  vector<int> centro_idx; // Centroider indices to handle
  std::deque<carma_obj<Tcomp> *> d_circularComs; //Circular buffer of commands for latency


  map<string, tuple<carma_obj<Tcomp> *, int, bool>> d_perturb_map;
  // perturbation command buffer

  carma_streams *streams;

  // allocation of d_centroids and d_com
  sutra_controller(carma_context *context, int nvalid, int nslope, int nactu,
                   float delay, sutra_dms *dms, int *idx_dms, int ndm,  int *idx_centro, int ncentro);
  virtual ~sutra_controller();

  virtual string get_type() = 0;

  //!!!! YOU MUST set d_centroids before calling it!!!!
  virtual int comp_com() = 0;

  // It is better to have something like this (+protected d_centroids):
  // virtual int comp_com (carma_obj<T> *new_centroids)=0;
  // it would imply copy, but would be much safer

  inline int nactu() { return this->nactus; }
  inline int nslope() { return this->nslopes; }

  cublasHandle_t cublas_handle() { return current_context->get_cublasHandle(); }

  void init_voltage() {
    init_voltage_impl<Tcomp, Tout>(this->d_voltage, this->d_comClipped);
  };

  int set_centroids_ref(Tcomp *centroids_ref);
  int add_perturb_voltage(string name, float *perturb, int N);
  int set_perturb_voltage(string name, float *perturb, int N);
  int remove_perturb_voltage(string name);
  int reset_perturb_voltage();
  int enable_perturb_voltage(string name);
  int disable_perturb_voltage(string name);
  int set_com(float *com, int nElem);
  int set_openloop(int open_loop_status, bool rst = true);
  int clip_commands();
  int comp_voltage();
  int comp_latency();
  int set_delay(float delay);
  int set_Vmin(float Vmin);
  int set_Vmax(float Vmax);
  int set_valMax(float valMax);
  int set_gain(float gain);
  int reset_coms();

  // int syevd_f(char meth, carma_obj<T> *d_U,
  //             carma_host_obj<T> *h_eingenvals);
  // int invgen(carma_obj<T> *d_mat, T cond, int job);
  int command_delay();
  int add_perturb();

 protected:
  mutex comp_voltage_mutex;
};

template <typename Tin, typename Tout>
typename std::enable_if<std::is_same<Tin, Tout>::value, void>::type
convertToVoltage(Tin *d_idata, Tout *d_odata, int N, float Vmin, float Vmax,
                 uint16_t valMax, carma_device *device){};

template <typename Tin, typename Tout>
typename std::enable_if<!std::is_same<Tin, Tout>::value, void>::type
convertToVoltage(Tin *d_idata, Tout *d_odata, int N, float Vmin, float Vmax,
                 uint16_t valMax, carma_device *device);

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
#endif  // _SUTRA_CONTROLLER_H_
