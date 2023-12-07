// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_controller_generic.h
//! \ingroup   libsutra
//! \class     SutraControllerGeneric
//! \brief     this class provides the controller_generic features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#ifndef _SUTRA_CONTROLLER_GENERIC_H_
#define _SUTRA_CONTROLLER_GENERIC_H_

#include <sutra_acquisim.h>
#include <sutra_controller.h>

template <typename Tcomp, typename Tout>
class SutraControllerGeneric : public SutraController<Tcomp, Tout> {
 public:
  CarmaObj<Tcomp> *d_matE;
  CarmaObj<Tcomp> *d_cmat;
  CarmaObj<Tcomp> *d_cmatPadded;
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
  int32_t fill_cmatPadded();
  int32_t distribute_cmat();

 private:
  template <typename Q = Tcomp>
  typename std::enable_if<!std::is_same<Q, half>::value, int32_t>::type
  fill_cmatPadded_impl() {
    return EXIT_SUCCESS;
  };
  template <typename Q = Tcomp>
  typename std::enable_if<std::is_same<Q, half>::value, int32_t>::type
  fill_cmatPadded_impl();
};

template <typename T>
void pad_cmat(T *idata, int32_t m, int32_t n, T *odata, int32_t m2, int32_t n2,
              CarmaDevice *device);
#endif  // _SUTRA_CONTROLLER_GENERIC_H_
