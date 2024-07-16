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

//! \file      sutra_controller.cpp
//! \ingroup   libsutra
//! \class     SutraController
//! \brief     this class provides the controller features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#include <sutra_controller.hpp>
#include <string>

template <typename Tcomp, typename Tout>
SutraController<Tcomp, Tout>::SutraController(CarmaContext *context,
                                                int32_t nslope, int32_t nactu, float delay,
                                                SutraDms *dms, int32_t *idx_dms,
                                                int32_t ndm, int32_t *idx_centro, int32_t ncentro) {
  this->current_context = context;
  this->device = context->get_active_device();

  this->nactus = nactu;
  this->nslopes = nslope;
  // current_context->set_active_device(device,1);
  this->d_com1 = nullptr;

  this->mainStream = NULL;
  cublasSetStream(this->cublas_handle(), this->mainStream);

  this->open_loop = 0;
  if (std::is_same<Tout, uint16_t>::value) {
    this->volt_min = -1.0f;
    this->volt_max = 1.0f;
  } else {
    this->volt_min = -1.e6;
    this->volt_max = 1.e6;
  }
  this->val_max = Tout(65535);

  this->set_delay(delay);

  int64_t dims_data1[2] = {1, 0};

  dims_data1[1] = nslope;
  this->d_centroids = new CarmaObj<Tcomp>(context, dims_data1);

  dims_data1[1] = nactu;
  this->d_com = new CarmaObj<Tcomp>(context, dims_data1);
  this->d_circular_coms.push_front(d_com);
  this->d_com1 = new CarmaObj<Tcomp>(context, dims_data1);
  if (this->delay > 0) {
    this->d_circular_coms.push_front(d_com1);
    while (this->d_circular_coms.size() <= int32_t(this->delay) + 1) {
      this->d_circular_coms.push_front(new CarmaObj<Tcomp>(context, dims_data1));
    }
  }
  this->d_com_clipped = new CarmaObj<Tcomp>(context, dims_data1);

  this->init_voltage();

  for (int32_t i = 0; i < ndm; i++) {
    this->d_dmseen.push_back(dms->d_dms[idx_dms[i]]);
  }

  for (int32_t i = 0; i < ncentro; i++) {
    this->centro_idx.push_back(idx_centro[i]);
  }

}

template <typename Tcomp, typename Tout>
int32_t SutraController<Tcomp, Tout>::set_gain(float gain) {
  this->gain = Tcomp(gain);
  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int32_t SutraController<Tcomp, Tout>::set_open_loop(int32_t open_loop_status,
                                                bool rst) {
  current_context->set_active_device(device, 1);
  this->open_loop = open_loop_status;

  if (this->open_loop && rst) {
    this->reset_coms();
  }
  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int32_t SutraController<Tcomp, Tout>::reset_coms() {
  current_context->set_active_device(device, 1);
  for (auto cobj : this->d_circular_coms) {
    cobj->reset();
  }

  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int32_t SutraController<Tcomp, Tout>::add_perturb_voltage(string name,
                                                       float *perturb, int32_t N) {
  std::lock_guard<std::mutex> lock(this->comp_voltage_mutex);
  current_context->set_active_device(device, 1);

  if (this->d_perturb_map.count(name) < 1) {
    int64_t dims_data2[3] = {2, N, this->nactu()};
    int32_t cpt = 0;
    CarmaObj<Tcomp> *d_perturb =
        new CarmaObj<Tcomp>(current_context, dims_data2);
    d_perturb->host2device(perturb);
    this->d_perturb_map[name] = std::make_tuple(d_perturb, cpt, true);
  } else
    DEBUG_TRACE("This perturb buffer already exists");

  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int32_t SutraController<Tcomp, Tout>::set_perturb_voltage(string name,
                                                       float *perturb, int32_t N) {
  current_context->set_active_device(device, 1);

  if (this->d_perturb_map.count(name)) {
    remove_perturb_voltage(name);
  }
  add_perturb_voltage(name, perturb, N);

  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int32_t SutraController<Tcomp, Tout>::remove_perturb_voltage(string name) {
  std::lock_guard<std::mutex> lock(this->comp_voltage_mutex);
  current_context->set_active_device(device, 1);
  typename map<string, tuple<CarmaObj<Tcomp> *, int32_t, bool>>::iterator it;

  if (this->d_perturb_map.count(name)) {
    it = this->d_perturb_map.find(name);
    delete std::get<0>(this->d_perturb_map[name]);
    this->d_perturb_map.erase(name);
  } else
    DEBUG_TRACE("This perturb buffer do not exist");

  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int32_t SutraController<Tcomp, Tout>::enable_perturb_voltage(string name) {
  std::lock_guard<std::mutex> lock(this->comp_voltage_mutex);
  current_context->set_active_device(device, 1);
  typename map<string, tuple<CarmaObj<Tcomp> *, int32_t, bool>>::iterator it;

  if (this->d_perturb_map.count(name)) {
    it = this->d_perturb_map.find(name);
    std::get<2>(this->d_perturb_map[name]) = true;
  } else
    DEBUG_TRACE("This perturb buffer do not exist");

  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int32_t SutraController<Tcomp, Tout>::disable_perturb_voltage(string name) {
  std::lock_guard<std::mutex> lock(this->comp_voltage_mutex);
  current_context->set_active_device(device, 1);
  typename map<string, tuple<CarmaObj<Tcomp> *, int32_t, bool>>::iterator it;

  if (this->d_perturb_map.count(name)) {
    it = this->d_perturb_map.find(name);
    std::get<2>(this->d_perturb_map[name]) = false;
  } else
    DEBUG_TRACE("This perturb buffer do not exist");

  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int32_t SutraController<Tcomp, Tout>::reset_perturb_voltage() {
  std::lock_guard<std::mutex> lock(this->comp_voltage_mutex);
  current_context->set_active_device(device, 1);
  typename map<string, tuple<CarmaObj<Tcomp> *, int32_t, bool>>::iterator it;
  it = this->d_perturb_map.begin();
  while (it != this->d_perturb_map.end()) {
    delete std::get<0>(it->second);
    it++;
  }
  this->d_perturb_map.clear();

  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int32_t SutraController<Tcomp, Tout>::command_delay() {
  current_context->set_active_device(device, 1);

  if (delay > 0) {
    CarmaObj<Tcomp> *tmp = std::move(this->d_circular_coms.front());
    tmp->copy(this->d_com, 1, 1);
    this->d_circular_coms.pop_front();
    this->d_circular_coms.push_back(tmp);
    this->d_com = this->d_circular_coms.back();
    this->d_com1 = this->d_circular_coms[this->d_circular_coms.size() - 2];
  }

  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int32_t SutraController<Tcomp, Tout>::clip_commands() {
  current_context->set_active_device(device, 1);
  if (this->volt_min < this->volt_max)
    this->d_com_clipped->clip(this->volt_min, this->volt_max, this->mainStream);
  else
    DEBUG_TRACE(
        "volt_max value must be greater than volt_min value. Nothing has been done.");
  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int32_t SutraController<Tcomp, Tout>::comp_latency() {
  // Command increment = a * d[k] + b*d[k-1] + c*d[k-2] + perturb
  this->current_context->set_active_device(this->device, 1);
  if (this->delay > 0) {
    this->d_com_clipped->reset();
    this->d_com_clipped->axpy(this->a, this->d_circular_coms[0], 1, 1);
    this->d_com_clipped->axpy(this->b, this->d_circular_coms[1], 1, 1);
  }

  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int32_t SutraController<Tcomp, Tout>::comp_voltage() {
  std::lock_guard<std::mutex> lock(this->comp_voltage_mutex);
  current_context->set_active_device(device, 1);
  if (this->delay > 0) {
    this->comp_latency();
    this->command_delay();
  } else {
    this->d_com_clipped->copy_from_async(this->d_com->get_data(),
                                 this->d_com_clipped->get_nb_elements(), this->mainStream);
  }
  this->add_perturb();
  this->clip_commands();
  convert_to_voltage<Tcomp, Tout>(
      this->d_com_clipped->get_data(), this->d_voltage->get_data(), this->nactu(),
      this->volt_min, this->volt_max, this->val_max,
      this->current_context->get_device(this->device), this->mainStream);

  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int32_t SutraController<Tcomp, Tout>::add_perturb() {
  typename map<string, tuple<CarmaObj<Tcomp> *, int32_t, bool>>::iterator it;
  int32_t cpt;
  CarmaObj<Tcomp> *d_perturb;
  int32_t any_perturb = 0;
  for (it = this->d_perturb_map.begin(); it != this->d_perturb_map.end();
       ++it) {
    if (std::get<2>(it->second)) {
      cpt = std::get<1>(it->second);
      d_perturb = std::get<0>(it->second);
      this->d_com_clipped->axpy(1.0f, d_perturb, d_perturb->get_dims(1), 1, cpt);
      if (cpt < d_perturb->get_dims(1) - 1)
        std::get<1>(it->second) = cpt + 1;
      else
        std::get<1>(it->second) = 0;
    }
  }

  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
SutraController<Tcomp, Tout>::~SutraController() {

  if (this->d_centroids != nullptr) delete this->d_centroids;
  for (int32_t k=0; k < this->d_circular_coms.size() ; k++) {
    delete this->d_circular_coms[k];
  }
  this->d_circular_coms.clear();
  if (this->d_com_clipped != nullptr &&
      (void *)this->d_com_clipped != (void *)this->d_com)
    delete this->d_com_clipped;
  if (this->d_voltage != nullptr &&
      (void *)this->d_voltage != (void *)this->d_com_clipped)
    delete this->d_voltage;
  this->reset_perturb_voltage();

}

template <typename Tcomp, typename Tout>
int32_t SutraController<Tcomp, Tout>::set_com(float *com, int32_t nElem) {
  current_context->set_active_device(device, 1);
  if (nElem == this->d_com->get_nb_elements()) {
    this->d_com->host2device(com);
    this->d_com_clipped->copy_from(this->d_com->get_data(),
                                 this->d_com_clipped->get_nb_elements());
  } else
    DEBUG_TRACE("Wrong dimension of com");

  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int32_t SutraController<Tcomp, Tout>::set_delay(float delay) {

  this->delay = delay;
  int32_t floor = (int32_t)delay;
  this->a = Tcomp(delay - floor);
  this->b = Tcomp(1 - a);


  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int32_t SutraController<Tcomp, Tout>::set_volt_min(float volt_min) {
  this->volt_min = volt_min;
  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int32_t SutraController<Tcomp, Tout>::set_volt_max(float volt_max) {
  this->volt_max = volt_max;
  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int32_t SutraController<Tcomp, Tout>::set_val_max(float val_max) {
  this->val_max = Tcomp(val_max);
  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int32_t SutraController<Tcomp, Tout>::comp_polc(CarmaObj<Tcomp>& sk,
   CarmaObj<Tcomp>& iMat, CarmaObj<Tcomp>& ol_meas){

  int64_t n_slope = iMat.get_dims(1);
  int64_t n_actu  = iMat.get_dims(2);

  std::string err = "";
  if(sk.get_dims(1) != n_slope){
    err += "\n\tinconsistent dimension for sk     : expected:";
    err += std::to_string(n_slope)+", got:" +std::to_string(sk.get_dims(1));
  }
  if(ol_meas.get_dims(1) != n_slope){
    err += "\n\tinconsistent dimension for ol_meas: expected:";
    err += std::to_string(n_slope)+", got:" +std::to_string(ol_meas.get_dims(1));
  }
  if(!err.empty()){
    throw std::runtime_error("ERROR while computing polc:" + err);
  }

  this->current_context->set_active_device(this->device, 1);
  // ol_meas = s_k
  ol_meas.copy(&sk, 1, 1);
  // ol_meas = - iMat * eff_u + s_k
  carma_gemv(this->cublas_handle(), 'n', n_slope, n_actu, Tcomp(-1.0f),
             iMat.get_data(), n_slope, d_com_clipped->get_data(), 1, Tcomp(1.0f),
             ol_meas.get_data(), 1);

  return EXIT_SUCCESS;
}

template class SutraController<float, float>;
template class SutraController<float, uint16_t>;
