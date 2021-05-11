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

//! \file      sutra_controller_generic_linear.cpp
//! \ingroup   libsutra
//! \class     sutra_controller_generic_linear
//! \brief     this class provides the controller_generic features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.1.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#include <sutra_controller_generic_linear.h>


template<typename T>
int rotate_circular_buffer(std::deque<CarmaObj<T> *> buffer){
  CarmaObj<T> *tmp = std::move(buffer.front());
    buffer.pop_front();
    buffer.push_back(tmp);
    return EXIT_SUCCESS;
}

template<typename T>
int update_circular(std::deque<CarmaObj<T> *> buffer, CarmaObj<T> &update){
    rotate_circular_buffer(buffer);
    buffer.at(0)->copy_from(update, update.get_nb_elements());
    return EXIT_SUCCESS;
}

template<typename T, typename Tout>
sutra_controller_generic_linear<T, Tout>::sutra_controller_generic_linear(
    CarmaContext *context, long nslopes, int nslope_buffers, long nactu, int nstates,
    int nstate_buffers, int nmodes, int nmode_buffers, int niir_in, int niir_out,
    float delay, bool polc, bool is_modal,
    SutraDms *dms, int *idx_dms, int ndm, int *idx_centro, int ncentro)
    : SutraController<T, Tout>(context, nslopes, nactu, delay, dms,
                                idx_dms, ndm, idx_centro, ncentro)
    {

  m_n_slope_buffers  = nslope_buffers;
  m_n_states  = nstates;
  m_n_state_buffers  = nstate_buffers;
  m_n_modes   = nmodes;
  m_n_mode_buffers   = nmode_buffers;
  m_n_iir_in  = niir_in;
  m_n_iir_out = niir_out;

  d_x_now = CarmaObj<T>(context, {1,m_n_states});
  for(int i=0;i<m_n_state_buffers;i++){
    d_circular_x.push_front(new CarmaObj<T>(context, {1,m_n_states}));
    d_A.push_back(CarmaObj<T>(context, {2,m_n_states, m_n_states}));
  }

  d_s_now = CarmaObj<T>(context, {1,nslopes});
  for(int i=0;i<m_n_slope_buffers;i++){
    d_circular_s.push_front(new CarmaObj<T>(context, {1,nslopes}));
    d_L.push_back(CarmaObj<T>(context, {2,m_n_states, nslopes}));
  }

  d_eff_u = CarmaObj<T>(context, {1,nactu});

  d_u_now = CarmaObj<T>(context, {1,nmodes});
  for(int i=0;i<m_n_iir_in;i++){
    d_circular_u_in.push_front(new CarmaObj<T>(context, {1,m_n_modes}));
    d_iir_b.push_back(CarmaObj<T>(context, {1,m_n_modes}));
  }
  for(int i=0;i<m_n_iir_out;i++){
    d_circular_u_out.push_front(new CarmaObj<T>(context, {1,m_n_modes}));
    d_iir_a.push_back(CarmaObj<T>(context, {1,m_n_modes}));
  }

  m_polc = false;
  m_modal = false;
  d_K = CarmaObj<T>(context, {2,nmodes , nstates});
  if(polc){
    m_polc = true;
    d_D = CarmaObj<T>(context, {2,nslopes , nactu  });
  }
  if(is_modal){
    m_modal = true;
    d_F = CarmaObj<T>(context, {2,nactu  , nmodes });
  }
  this->current_context->set_active_device(this->device, 1);
  //cudaEventCreateWithFlags(&start_mvm_event, cudaEventDisableTiming);
}

template<typename T, typename Tout>
sutra_controller_generic_linear<T, Tout>::~sutra_controller_generic_linear() {
  this->current_context->set_active_device(this->device, 1);
  for(auto &c : d_circular_x)    {delete c;}
  for(auto &c : d_circular_s)    {delete c;}
  for(auto &c : d_circular_u_in) {delete c;}
  for(auto &c : d_circular_u_out){delete c;}
  d_circular_x.clear();
  d_circular_s.clear();
  d_circular_u_in.clear();
  d_circular_u_out.clear();
}



template <typename T, typename Tout>
string sutra_controller_generic_linear<T, Tout>::get_type() {
  this->current_context->set_active_device(this->device, 1);
  return "generic_linear";
}

template<typename T, typename Tout>
int sutra_controller_generic_linear<T, Tout>::set_polc(bool p) {
  m_polc = p;
  return EXIT_SUCCESS;
}

template<typename T, typename Tout>
int sutra_controller_generic_linear<T, Tout>::set_A(float *M, int i) {
  this->current_context->set_active_device(this->device, 1);
  d_A[i].host2device(M);
  return EXIT_SUCCESS;
}

template<typename T, typename Tout>
int sutra_controller_generic_linear<T, Tout>::set_L(float *M, int i) {
  this->current_context->set_active_device(this->device, 1);
  d_L[i].host2device(M);
  return EXIT_SUCCESS;
}

template<typename T, typename Tout>
int sutra_controller_generic_linear<T, Tout>::set_K(float *M) {
  this->current_context->set_active_device(this->device, 1);
  d_K.host2device(M);
  return EXIT_SUCCESS;
}

template<typename T, typename Tout>
int sutra_controller_generic_linear<T, Tout>::set_D(float *M) {
  this->current_context->set_active_device(this->device, 1);
  d_D.host2device(M);
  return EXIT_SUCCESS;
}

template<typename T, typename Tout>
int sutra_controller_generic_linear<T, Tout>::set_F(float *M) {
  this->current_context->set_active_device(this->device, 1);
  d_F.host2device(M);
  return EXIT_SUCCESS;
}

template<typename T, typename Tout>
int sutra_controller_generic_linear<T, Tout>::set_iir_a(float *M, int i) {
  this->current_context->set_active_device(this->device, 1);
  d_iir_a[i].host2device(M);
  return EXIT_SUCCESS;
}

template<typename T, typename Tout>
int sutra_controller_generic_linear<T, Tout>::set_iir_b(float *M, int i) {
  this->current_context->set_active_device(this->device, 1);
  d_iir_b[i].host2device(M);
  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int sutra_controller_generic_linear<T, Tout>::comp_polc(){
  comp_polc(a, *(d_circular_coms.at(1)), b, *(d_circular_coms.at(0)),
    *(d_centroids), d_D, d_s_now, d_eff_u);

  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int sutra_controller_generic_linear<T, Tout>::comp_com() {
  this->current_context->set_active_device(this->device, 1);
  if(m_polc){
    comp_polc();
  }else{
    d_s_now.copy_from(*d_centroids, d_centroids->get_nb_elements());
  }
  update_circular(d_circular_s, d_s_now);

  recursion();
  innovation();
  modal_projection();

  update_circular(d_circular_x, d_x_now);

  update_circular(d_circular_u_in, d_u_now);

  filter_iir_in();
  filter_iir_out();

  update_circular(d_circular_u_out, d_u_now);

  if(m_modal){
    carma_gemv(cublas_handle(), 'n', nactus, m_n_modes, T(1.0f),
               d_F.get_data(), m_n_modes, d_u_now.get_data(), 1, T(1.0f),
               d_com->get_data(), 1);
  }else{
    d_com->copy_from(d_u_now, d_u_now.get_nb_elements());
  }

  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int sutra_controller_generic_linear<T, Tout>::recursion(){
  //x_now = sum_{i}^{m_n_state_buffer}(A[i]* circular_x[i])
  if(d_circular_x.size()>0){
    //1st gemv : reset x_now
    carma_gemv(cublas_handle(), 'n', m_n_states, m_n_states, T(1.0f),
               d_A[0].get_data(), m_n_states, d_circular_x.at(0)->get_data(), 1, T(0.0f),
               d_x_now.get_data(), 1);
    //continue summation
    for(int i=1;i<d_circular_x.size(); i++){
      carma_gemv(cublas_handle(), 'n', m_n_states, m_n_states, T(1.0f),
                 d_A[i].get_data(), m_n_states, d_circular_x.at(i)->get_data(), 1, T(1.0f),
                 d_x_now.get_data(), 1);
    }
  }
  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int sutra_controller_generic_linear<T, Tout>::innovation(){
  //x_now += sum_{i}^{m_n_slope_buffer}(L[i]* circular_s[i])
  for(int i=0;i<d_circular_s.size(); i++){
    carma_gemv(cublas_handle(), 'n', m_n_states, nslopes, T(1.0f),
               d_L[i].get_data(), m_n_states, d_circular_s.at(i)->get_data(), 1, T(1.0f),
               d_x_now.get_data(), 1);
  }
  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int sutra_controller_generic_linear<T, Tout>::modal_projection(){
  if(m_n_state_buffers>0){
    // u_now = K * x_now
    carma_gemv(cublas_handle(), 'n', m_n_modes, m_n_states, T(1.0f),
               d_K.get_data(), m_n_modes, d_x_now.get_data(), 1, T(1.0f),
               d_u_now.get_data(), 1);
  }else{
    //u_now = x_now
    d_u_now.copy_from(d_x_now,d_x_now.get_nb_elements());
  }
  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int sutra_controller_generic_linear<T, Tout>::filter_iir_in(){
   // u_now = sum_{i=0}^{}(iir_b[i]*u_in[i]) //gbmv
  if(m_n_iir_in>0){
    //1st gemv : reset x_now
    carma_sbmv(cublas_handle(), 'l', m_n_modes, 0, T(1.0f),
               d_iir_b[0].get_data(), 1, d_circular_u_in.at(0)->get_data(), 1, T(0.0f),
               d_u_now.get_data(), 1);

    //continue summation
    for(int i=1;i<m_n_iir_in; i++){
    carma_sbmv(cublas_handle(), 'l', m_n_modes, 0, T(1.0f),
               d_iir_b[i].get_data(), 1, d_circular_u_in.at(i)->get_data(), 1, T(1.0f),
               d_u_now.get_data(), 1);
    }
  }
  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int sutra_controller_generic_linear<T, Tout>::filter_iir_out(){
   // u_now += sum_{i=0}^{}(iir_a[i]*u_out[i]) //gbmv
    for(int i=1;i<m_n_iir_in; i++){
    carma_sbmv(cublas_handle(), 'l', m_n_modes, 0, T(1.0f),
               d_iir_a[i].get_data(), 1, d_circular_u_out.at(i)->get_data(), 1, T(1.0f),
               d_u_now.get_data(), 1);
  }

  return EXIT_SUCCESS;
}


template class sutra_controller_generic_linear<float, float>;
template class sutra_controller_generic_linear<float, uint16_t>;
