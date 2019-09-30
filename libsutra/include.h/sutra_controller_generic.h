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

//! \file      sutra_controller_generic.h
//! \ingroup   libsutra
//! \class     sutra_controller_generic
//! \brief     this class provides the controller_generic features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.3.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#ifndef _SUTRA_CONTROLLER_GENERIC_H_
#define _SUTRA_CONTROLLER_GENERIC_H_

#include <sutra_acquisim.h>
#include <sutra_controller.h>

template <typename Tcomp, typename Tout>
class sutra_controller_generic : public sutra_controller<Tcomp, Tout> {
 public:
  carma_obj<Tcomp> *d_matE;
  carma_obj<Tcomp> *d_cmat;
  carma_obj<Tcomp> *d_cmatPadded;
  carma_obj<Tcomp> *d_gain;
  carma_obj<Tcomp> *d_decayFactor;
  carma_obj<Tcomp> *d_compbuff;  // Buffer for computations
  carma_obj<Tcomp> *d_compbuff2;
  carma_obj<Tcomp> *d_olmeas;  // Open-loop measurements for POLC
  carma_obj<Tcomp> *d_imat;

  bool polc;
  int nstates;

  string command_law;

 public:
  sutra_controller_generic(carma_context *context, long nvalid, long nslope,
                           long nactu, float delay, sutra_dms *dms,
                           int *idx_dms, int ndm, int *idx_centro, int ncentro, int nstates);
  sutra_controller_generic(const sutra_controller_generic &controller);
  ~sutra_controller_generic();

  string get_type();
  string get_commandlaw();
  int set_decayFactor(float *decayFactor);
  int set_mgain(float *gain);
  int set_cmat(float *cmat);
  int set_matE(float *matE);
  int set_commandlaw(string law);
  int set_polc(bool p);
  int set_imat(float *imat);
  int comp_polc();
  int comp_com();
  int fill_cmatPadded();

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
              carma_device *device);
#endif  // _SUTRA_CONTROLLER_GENERIC_H_
