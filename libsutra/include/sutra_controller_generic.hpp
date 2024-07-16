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

//! \file      sutra_controller_generic.hpp
//! \ingroup   libsutra
//! \class     SutraControllerGeneric
//! \brief     this class provides the controller_generic features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#ifndef _SUTRA_CONTROLLER_GENERIC_H_
#define _SUTRA_CONTROLLER_GENERIC_H_

#include <sutra_acquisim.hpp>
#include <sutra_controller.hpp>

template <typename Tcomp, typename Tout>
class SutraControllerGeneric : public SutraController<Tcomp, Tout> {
 public:
  CarmaObj<Tcomp> *d_matE;
  CarmaObj<Tcomp> *d_cmat;
  CarmaObj<Tcomp> *d_gain;
  CarmaObj<Tcomp> *d_decayFactor;
  CarmaObj<Tcomp> *d_compbuff;  // Buffer for computations
  CarmaObj<Tcomp> *d_compbuff2;
  CarmaObj<Tcomp> *d_olmeas;  // Open-loop measurements for POLC
  CarmaObj<Tcomp> *d_imat;
  std::vector<CarmaObj<Tcomp> *> d_err_ngpu;
  std::vector<CarmaObj<Tcomp> *> d_centroids_ngpu;
  std::vector<CarmaObj<Tcomp> *> d_cmat_ngpu;
  std::vector<int32_t> P2Pdevices;
  std::vector<cudaEvent_t> events;
  cudaEvent_t start_mvm_event;
  std::vector<cudaStream_t> streams;
  bool polc;
  int32_t nstates;
  Tcomp leaky_factor;
  string command_law;

 public:
  SutraControllerGeneric(CarmaContext *context, int64_t nslope,
                           int64_t nactu, float delay, SutraDms *dms,
                           int32_t *idx_dms, int32_t ndm, int32_t *idx_centro, int32_t ncentro, int32_t nstates);
  SutraControllerGeneric(const SutraControllerGeneric &controller);
  ~SutraControllerGeneric();

  string get_type();
  string get_commandlaw();
  int32_t set_decayFactor(float *decayFactor);
  int32_t set_modal_gains(float *gain);
  int32_t set_cmat(float *cmat);
  int32_t set_matE(float *matE);
  int32_t set_commandlaw(string law);
  int32_t set_polc(bool p);
  int32_t set_imat(float *imat);
  int32_t set_leaky_factor(Tcomp factor);
  using SutraController<Tcomp,Tout>::comp_polc;
  int32_t comp_polc();
  int32_t comp_com();
  int32_t distribute_cmat();
};

#endif  // _SUTRA_CONTROLLER_GENERIC_H_
