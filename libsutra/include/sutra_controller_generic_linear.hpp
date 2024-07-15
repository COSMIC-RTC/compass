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

//! \file      sutra_controller_generic_linear.hpp
//! \ingroup   libsutra
//! \class     SutraControllerGenericLinear
//! \brief     this class provides the controller_generic features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \version   5.5.0
//! \date      2022/01/24

#ifndef _SUTRA_CONTROLLER_GENERIC_LINEAR_H_
#define _SUTRA_CONTROLLER_GENERIC_LINEAR_H_

#include <sutra_acquisim.hpp>
#include <sutra_controller.hpp>

#include <memory>

template <typename Tcomp, typename Tout>
class SutraControllerGenericLinear : public SutraController<Tcomp, Tout> {


private:
  using SutraController<Tcomp,Tout>::nslopes; //!< num of meas in slope vector
  using SutraController<Tcomp,Tout>::nactus;  //!< num of coms in com vector

  bool m_polc;              //!< flag to compute POL slopes
  bool m_modal;             //!< flag to do projection from modes to actu
  int32_t m_n_slope_buffers;    //!< num of historic slopes vectors to use
  int32_t m_n_states;           //!< num of states in state vector
  int32_t m_n_state_buffers;    //!< num of historic state vectors to use
  int32_t m_n_modes;            //!< num of modes in mode vector
  int32_t m_n_iir_in;           //!< num of input mode vectors for iir filter
  int32_t m_n_iir_out;          //!< num of output mode vectors for iir filter

  using SutraController<Tcomp,Tout>::a; //!< coefficient used to comp delay
  using SutraController<Tcomp,Tout>::b; //!< coefficient used to comp delay

  using SutraController<Tcomp,Tout>::cublas_handle;

public:

  bool polc(){return m_polc;}
  bool modal(){return m_modal;}
  int32_t n_slope_buffers(){return m_n_slope_buffers;}
  int32_t n_states(){return m_n_states;}
  int32_t n_state_buffers(){return m_n_state_buffers;}
  int32_t n_modes(){return m_n_modes;}
  int32_t n_iir_in(){return m_n_iir_in;}
  int32_t n_iir_out(){return m_n_iir_out;}

  std::unique_ptr<CarmaObj<Tcomp>> d_x_now;   //!< vector used for state calcs
  std::unique_ptr<CarmaObj<Tcomp>> d_s_now;   //!< vector used for slope calcs
  std::unique_ptr<CarmaObj<Tcomp>> d_u_now;   //!< vector used for modal calcs

  std::deque<CarmaObj<Tcomp> *> d_circular_x;     //!< circ buffer states
  std::deque<CarmaObj<Tcomp> *> d_circular_s;     //!< circ buffer slopes
  std::deque<CarmaObj<Tcomp> *> d_circular_u_in;  //!< circ buffer iir inputs
  std::deque<CarmaObj<Tcomp> *> d_circular_u_out; //!< circ buffer iir output

  std::vector<CarmaObj<Tcomp>*> d_matA;      //!< list of A matrices (recursions)
  std::vector<CarmaObj<Tcomp>*> d_matL;      //!< list of L matrices (innovations)
  std::unique_ptr<CarmaObj<Tcomp>> d_matK;   //!< K matrix (state to modes)
  std::unique_ptr<CarmaObj<Tcomp>> d_matD;   //!< D matrix (interaction matrix)
  std::unique_ptr<CarmaObj<Tcomp>> d_matF;   //!< F matrix (mode to actu)

  std::vector<CarmaObj<Tcomp>*> d_iir_a;  //!< list of iir 'a' vectors (outputs)
  std::vector<CarmaObj<Tcomp>*> d_iir_b;  //!< list of iir 'b' vectors (inputs)

  using SutraController<Tcomp,Tout>::d_com;  //!< most recently computed command
  using SutraController<Tcomp,Tout>::d_centroids;  //!< closed loop slope vector
  using SutraController<Tcomp,Tout>::d_circular_coms; //!< circular com buffer
  using SutraController<Tcomp,Tout>::delay;

 public:
 SutraControllerGenericLinear() = delete;
  SutraControllerGenericLinear(
    CarmaContext *context, int32_t nslope, int32_t nslopes_buffers, int32_t nactu, int32_t nstates,
    int32_t nstates_buffers, int32_t nmodes, int32_t niir_in, int32_t niir_out,
    float delay, bool polc, bool is_modal,
    SutraDms *dms, int32_t *idx_dms, int32_t ndm, int32_t *idx_centro, int32_t ncentro);
  //SutraControllerGenericLinear(const SutraControllerGenericLinear &controller);
  ~SutraControllerGenericLinear();

  string get_type();
  string get_commandlaw();
  int32_t set_matA(float * A, int32_t i);
  int32_t set_matL(float * L, int32_t i);
  int32_t set_matK(float * K);
  int32_t set_matD(float * D);
  int32_t set_matF(float * F);
  int32_t set_iir_a(float * iir_a, int32_t i);
  int32_t set_iir_b(float * iir_b, int32_t i);
  int32_t set_polc(bool p);
  using SutraController<Tcomp,Tout>::comp_polc;
  int32_t comp_polc();
  int32_t comp_com();

  int32_t recursion();
  int32_t innovation();
  int32_t modal_projection();
  int32_t filter_iir_in();
  int32_t filter_iir_out();

};

#endif  // _SUTRA_CONTROLLER_GENERIC_LINEAR_H_
