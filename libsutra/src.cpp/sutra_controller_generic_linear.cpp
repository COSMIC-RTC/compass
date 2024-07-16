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

//! \file      sutra_controller_generic_linear.cpp
//! \ingroup   libsutra
//! \class     SutraControllerGenericLinear
//! \brief     this class provides the controller_generic features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#include <sutra_controller_generic_linear.hpp>

template<typename T>
int32_t rotate_circular_buffer(std::deque<CarmaObj<T> *> &buffer){
  CarmaObj<T> *tmp = std::move(buffer.back());
  buffer.pop_back();
  buffer.push_front(tmp);
  return EXIT_SUCCESS;
}

template<typename T>
int32_t update_circular(std::deque<CarmaObj<T> *> &buffer, CarmaObj<T> &update){
  if(buffer.size()>0){
    rotate_circular_buffer(buffer);
    buffer.front()->copy_from(update, update.get_nb_elements());
  }
    return EXIT_SUCCESS;
}

template<typename T>
std::unique_ptr<CarmaObj<T>> make_unique_carma_obj(CarmaContext *context, std::initializer_list<int64_t> il) {

  std::vector<int64_t> dims(il.size() + 1);
  dims[0] = il.size();

  std::copy(il.begin(), il.end(), dims.begin() + 1);

  return std::make_unique<CarmaObj<T>>(context, dims);
}

template<typename T, typename Tout>
SutraControllerGenericLinear<T, Tout>::SutraControllerGenericLinear(
    CarmaContext *context, int32_t nslope, int32_t nslope_buffers, int32_t nactu, int32_t nstates,
    int32_t nstate_buffers, int32_t nmodes, int32_t niir_in, int32_t niir_out,
    float delay, bool polc, bool is_modal,
    SutraDms *dms, int32_t *idx_dms, int32_t ndm, int32_t *idx_centro, int32_t ncentro)
    : SutraController<T, Tout>(context, nslope, nactu, delay, dms,
                                idx_dms, ndm, idx_centro, ncentro){
  // TODO: remove n_mode_buffers everywhere
  m_n_slope_buffers  = nslope_buffers;
  m_n_states  = nstates;
  m_n_state_buffers  = nstate_buffers;
  m_n_modes   = nmodes;
  m_n_iir_in  = niir_in;
  m_n_iir_out = niir_out;

  d_x_now = make_unique_carma_obj<T>(context, {m_n_states});
  for(int32_t i=0;i<m_n_state_buffers;i++){
    d_circular_x.push_front(new CarmaObj<T>(context, {1,m_n_states}));
    d_matA.push_back(new CarmaObj<T>(context, {2,m_n_states, m_n_states}));
  }
  d_s_now = make_unique_carma_obj<T>(context, {nslopes});
  for(int32_t i=0;i<m_n_slope_buffers;i++){
    d_circular_s.push_front(new CarmaObj<T>(context, {1,nslopes}));
    d_matL.push_back(new CarmaObj<T>(context, {2,m_n_states, nslopes}));
  }

  d_u_now = make_unique_carma_obj<T>(context, {nmodes});
  for(int32_t i=0;i<m_n_iir_in;i++){
    d_circular_u_in.push_front(new CarmaObj<T>(context, {1,m_n_modes}));
    d_iir_b.push_back(new CarmaObj<T>(context, {1,m_n_modes}));
  }
  for(int32_t i=0;i<m_n_iir_out;i++){
    d_circular_u_out.push_front(new CarmaObj<T>(context, {1,m_n_modes}));
    d_iir_a.push_back(new CarmaObj<T>(context, {1,m_n_modes}));
  }

  m_polc = false;
  m_modal = false;
  d_matK = make_unique_carma_obj<T>(context, {nmodes , nstates});
  if(polc){
    m_polc = true;
    d_matD = make_unique_carma_obj<T>(context, {nslopes , nactu  });
  }
  if(is_modal){
    m_modal = true;
    d_matF = make_unique_carma_obj<T>(context, {nactu  , nmodes });
  }
  this->current_context->set_active_device(this->device, 1);
}

template<typename T, typename Tout>
SutraControllerGenericLinear<T, Tout>::~SutraControllerGenericLinear() {
  this->current_context->set_active_device(this->device, 1);
  for(auto &c : d_circular_x)    {delete c;}
  for(auto &c : d_circular_s)    {delete c;}
  for(auto &c : d_circular_u_in) {delete c;}
  for(auto &c : d_circular_u_out){delete c;}
  for(auto &v : d_iir_a){delete v;}
  for(auto &v : d_iir_b){delete v;}
  for(auto &v : d_matA){delete v;}
  for(auto &v : d_matL){delete v;}
  d_circular_x.clear();
  d_circular_s.clear();
  d_circular_u_in.clear();
  d_circular_u_out.clear();
}

template <typename T, typename Tout>
string SutraControllerGenericLinear<T, Tout>::get_type() {
  this->current_context->set_active_device(this->device, 1);
  return "generic_linear";
}

template<typename T, typename Tout>
int32_t SutraControllerGenericLinear<T, Tout>::set_polc(bool p) {
  m_polc = p;
  return EXIT_SUCCESS;
}

template<typename T, typename Tout>
int32_t SutraControllerGenericLinear<T, Tout>::set_matA(float *M, int32_t i) {
  this->current_context->set_active_device(this->device, 1);
  d_matA[i]->host2device(M);
  return EXIT_SUCCESS;
}

template<typename T, typename Tout>
int32_t SutraControllerGenericLinear<T, Tout>::set_matL(float *M, int32_t i) {
  this->current_context->set_active_device(this->device, 1);
  d_matL[i]->host2device(M);
  return EXIT_SUCCESS;
}

template<typename T, typename Tout>
int32_t SutraControllerGenericLinear<T, Tout>::set_matK(float *M) {
  this->current_context->set_active_device(this->device, 1);
  (*d_matK).host2device(M);
  return EXIT_SUCCESS;
}

template<typename T, typename Tout>
int32_t SutraControllerGenericLinear<T, Tout>::set_matD(float *M) {
  this->current_context->set_active_device(this->device, 1);
  (*d_matD).host2device(M);
  return EXIT_SUCCESS;
}

template<typename T, typename Tout>
int32_t SutraControllerGenericLinear<T, Tout>::set_matF(float *M) {
  this->current_context->set_active_device(this->device, 1);
  (*d_matF).host2device(M);
  return EXIT_SUCCESS;
}

template<typename T, typename Tout>
int32_t SutraControllerGenericLinear<T, Tout>::set_iir_a(float *M, int32_t i) {
  this->current_context->set_active_device(this->device, 1);
  d_iir_a[i]->host2device(M);
  return EXIT_SUCCESS;
}

template<typename T, typename Tout>
int32_t SutraControllerGenericLinear<T, Tout>::set_iir_b(float *M, int32_t i) {
  this->current_context->set_active_device(this->device, 1);
  d_iir_b[i]->host2device(M);
  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int32_t SutraControllerGenericLinear<T, Tout>::comp_polc(){
  comp_polc(*d_centroids, *d_matD, *d_s_now);
  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int32_t SutraControllerGenericLinear<T, Tout>::comp_com() {
  this->current_context->set_active_device(this->device, 1);
  if(m_polc){
    comp_polc();
  }else{
    (*d_s_now).copy_from(*d_centroids, d_centroids->get_nb_elements());
  }
  update_circular(d_circular_s, *d_s_now);

  recursion();
  innovation();
  modal_projection();

  update_circular(d_circular_x, *d_x_now);

  update_circular(d_circular_u_in, *d_u_now);

  filter_iir_in();
  filter_iir_out();

  update_circular(d_circular_u_out, *d_u_now);

  if(m_modal){
    carma_gemv(cublas_handle(), 'n', nactus, m_n_modes, T(1.0f),
               (*d_matF).get_data(), nactus, (*d_u_now).get_data(), 1, T(0.0f),
               d_com->get_data(), 1);
  }else{
    d_com->copy_from(*d_u_now, (*d_u_now).get_nb_elements());
  }

  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int32_t SutraControllerGenericLinear<T, Tout>::recursion(){
  //x_now = sum_{i}^{m_n_state_buffer}(A[i]* circular_x[i])
  if(d_circular_x.size()>0){
    //1st gemv : reset x_now
    carma_gemv(cublas_handle(), 'n', m_n_states, m_n_states, T(1.0f),
               d_matA[0]->get_data(), m_n_states, d_circular_x.at(0)->get_data(), 1, T(0.0f),
               (*d_x_now).get_data(), 1);
    //continue summation
    for(int32_t i=1;i<d_circular_x.size(); i++){
      carma_gemv(cublas_handle(), 'n', m_n_states, m_n_states, T(1.0f),
                 d_matA[i]->get_data(), m_n_states, d_circular_x.at(i)->get_data(), 1, T(1.0f),
                 (*d_x_now).get_data(), 1);
    }
  } else {
    (*d_x_now).memset(T(0.0f));
  }
  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int32_t SutraControllerGenericLinear<T, Tout>::innovation(){
  //x_now += sum_{i}^{m_n_slope_buffer}(L[i]* circular_s[i])
  for(int32_t i=0;i<d_circular_s.size(); i++){
    carma_gemv(cublas_handle(), 'n', m_n_states, nslopes, T(1.0f),
               d_matL[i]->get_data(), m_n_states, d_circular_s.at(i)->get_data(), 1, T(1.0f),
               (*d_x_now).get_data(), 1);
  }
  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int32_t SutraControllerGenericLinear<T, Tout>::modal_projection(){
  if(m_n_state_buffers>0){
    // u_now = K * x_now
    carma_gemv(cublas_handle(), 'n', m_n_modes, m_n_states, T(1.0f),
               (*d_matK).get_data(), m_n_modes, (*d_x_now).get_data(), 1, T(0.0f),
               (*d_u_now).get_data(), 1);
  }else{
    // u_now = x_now
    (*d_u_now).copy_from(*d_x_now,(*d_x_now).get_nb_elements());
  }
  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int32_t SutraControllerGenericLinear<T, Tout>::filter_iir_in(){
   // u_now = sum_{i=0}^{}(iir_b[i]*u_in[i]) //gbmv
  if(m_n_iir_in>0){
    //1st gemv : reset u_now
    carma_sbmv(cublas_handle(), 'l', m_n_modes, 0, T(1.0f),
               d_iir_b[0]->get_data(), 1, d_circular_u_in.at(0)->get_data(), 1, T(0.0f),
               (*d_u_now).get_data(), 1);

    //continue summation
    for(int32_t i=1;i<m_n_iir_in; i++){
    carma_sbmv(cublas_handle(), 'l', m_n_modes, 0, T(1.0f),
               d_iir_b[i]->get_data(), 1, d_circular_u_in.at(i)->get_data(), 1, T(1.0f),
               (*d_u_now).get_data(), 1);
    }
  }
  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int32_t SutraControllerGenericLinear<T, Tout>::filter_iir_out(){
   // u_now += sum_{i=0}^{}(iir_a[i]*u_out[i]) //gbmv
  for(int32_t i=0;i<m_n_iir_out; i++){
    carma_sbmv(cublas_handle(), 'l', m_n_modes, 0, T(1.0f),
               d_iir_a[i]->get_data(), 1, d_circular_u_out.at(i)->get_data(), 1, T(1.0f),
               (*d_u_now).get_data(), 1);
  }
  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int32_t SutraController<Tcomp, Tout>::reset_coms() {
  current_context->set_active_device(device, 1);
  for (auto cobj : this->d_circular_coms) {
    cobj->reset();
  }
  for (auto cobj : this->d_circular_x) {
    cobj->reset();
  }
  for (auto cobj : this->d_circular_s) {
    cobj->reset();
  }
  for (auto cobj : this->d_circular_u_in) {
    cobj->reset();
  }
  for (auto cobj : this->d_circular_u_out) {
    cobj->reset();
  }
  return EXIT_SUCCESS;
}

template class SutraControllerGenericLinear<float, float>;
template class SutraControllerGenericLinear<float, uint16_t>;
