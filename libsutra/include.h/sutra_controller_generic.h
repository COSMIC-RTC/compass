// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutra_controller_generic.h
//! \ingroup   libsutra
//! \class     sutra_controller_generic
//! \brief     this class provides the controller_generic features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.3.0
//! \date      2022/01/24

#ifndef _SUTRA_CONTROLLER_GENERIC_H_
#define _SUTRA_CONTROLLER_GENERIC_H_

#include <sutra_acquisim.h>
#include <sutra_controller.h>

template <typename Tcomp, typename Tout>
class sutra_controller_generic : public SutraController<Tcomp, Tout> {
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
  std::vector<int> P2Pdevices;
  std::vector<cudaEvent_t> events;
  cudaEvent_t start_mvm_event;
  std::vector<cudaStream_t> streams;
  bool polc;
  int nstates;
  Tcomp leaky_factor;
  string command_law;

 public:
  sutra_controller_generic(CarmaContext *context, long nslope,
                           long nactu, float delay, SutraDms *dms,
                           int *idx_dms, int ndm, int *idx_centro, int ncentro, int nstates);
  sutra_controller_generic(const sutra_controller_generic &controller);
  ~sutra_controller_generic();

  string get_type();
  string get_commandlaw();
  int set_decayFactor(float *decayFactor);
  int set_modal_gains(float *gain);
  int set_cmat(float *cmat);
  int set_matE(float *matE);
  int set_commandlaw(string law);
  int set_polc(bool p);
  int set_imat(float *imat);
  int set_leaky_factor(Tcomp factor);
  using SutraController<Tcomp,Tout>::comp_polc;
  int comp_polc();
  int comp_com();
  int fill_cmatPadded();
  int distribute_cmat();

 private:
  template <typename Q = Tcomp>
  typename std::enable_if<!std::is_same<Q, half>::value, int>::type
  fill_cmatPadded_impl() {
    return EXIT_SUCCESS;
  };
  template <typename Q = Tcomp>
  typename std::enable_if<std::is_same<Q, half>::value, int>::type
  fill_cmatPadded_impl();
};

template <typename T>
void pad_cmat(T *idata, int m, int n, T *odata, int m2, int n2,
              CarmaDevice *device);
#endif  // _SUTRA_CONTROLLER_GENERIC_H_
