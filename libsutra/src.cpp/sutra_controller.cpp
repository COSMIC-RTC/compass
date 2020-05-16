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

//! \file      sutra_controller.cpp
//! \ingroup   libsutra
//! \class     SutraController
//! \brief     this class provides the controller features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.0.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#include <sutra_controller.h>
#include <string>

template <typename Tcomp, typename Tout>
SutraController<Tcomp, Tout>::SutraController(CarmaContext *context,
                                                int nvalid, int nslope,
                                                int nactu, float delay,
                                                SutraDms *dms, int *idx_dms,
                                                int ndm, int *idx_centro, int ncentro) {
  this->current_context = context;
  this->device = context->get_active_device();

  this->nactus = nactu;
  this->nslopes = nslope;
  // current_context->set_active_device(device,1);
  this->d_com1 = nullptr;
  this->d_com_padded = nullptr;
  this->d_centroids_padded = nullptr;

  int nstreams = 1;  // nvalid/10;

  while (nactu % nstreams != 0) nstreams--;

  std::cerr << "controller uses " << nstreams << " streams" << std::endl;
  streams = new CarmaStreams(nstreams);

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

  long dims_data1[2] = {1, 0};

  dims_data1[1] = nslope;
  this->d_centroids = new CarmaObj<Tcomp>(context, dims_data1);

  dims_data1[1] = nactu;
  this->d_com = new CarmaObj<Tcomp>(context, dims_data1);
  this->d_circular_coms.push_front(d_com);
  this->d_com1 = new CarmaObj<Tcomp>(context, dims_data1);
  if (this->delay > 0) {
    this->d_circular_coms.push_front(d_com1);
    while (this->d_circular_coms.size() <= int(this->delay) + 1) {
      this->d_circular_coms.push_front(new CarmaObj<Tcomp>(context, dims_data1));
    }
  }
  this->d_com_clipped = new CarmaObj<Tcomp>(context, dims_data1);

  this->init_voltage();

  for (int i = 0; i < ndm; i++) {
    this->d_dmseen.push_back(dms->d_dms[idx_dms[i]]);
  }

  for (int i = 0; i < ncentro; i++) {
    this->centro_idx.push_back(idx_centro[i]);
  }

}

template <typename Tcomp, typename Tout>
int SutraController<Tcomp, Tout>::set_gain(float gain) {
  this->gain = Tcomp(gain);
  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int SutraController<Tcomp, Tout>::set_open_loop(int open_loop_status,
                                                bool rst) {
  current_context->set_active_device(device, 1);
  this->open_loop = open_loop_status;

  if (this->open_loop && rst) {
    this->reset_coms();
  }
  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int SutraController<Tcomp, Tout>::reset_coms() {
  current_context->set_active_device(device, 1);
  for (auto cobj : this->d_circular_coms) {
    cobj->reset();
  }

  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int SutraController<Tcomp, Tout>::add_perturb_voltage(string name,
                                                       float *perturb, int N) {
  std::lock_guard<std::mutex> lock(this->comp_voltage_mutex);
  current_context->set_active_device(device, 1);

  if (this->d_perturb_map.count(name) < 1) {
    long dims_data2[3] = {2, N, this->nactu()};
    int cpt = 0;
    CarmaObj<Tcomp> *d_perturb =
        new CarmaObj<Tcomp>(current_context, dims_data2);
    d_perturb->host2device(perturb);
    this->d_perturb_map[name] = std::make_tuple(d_perturb, cpt, true);
  } else
    DEBUG_TRACE("This perturb buffer already exists");

  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int SutraController<Tcomp, Tout>::set_perturb_voltage(string name,
                                                       float *perturb, int N) {
  current_context->set_active_device(device, 1);

  if (this->d_perturb_map.count(name)) {
    remove_perturb_voltage(name);
  }
  add_perturb_voltage(name, perturb, N);

  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int SutraController<Tcomp, Tout>::remove_perturb_voltage(string name) {
  std::lock_guard<std::mutex> lock(this->comp_voltage_mutex);
  current_context->set_active_device(device, 1);
  typename map<string, tuple<CarmaObj<Tcomp> *, int, bool>>::iterator it;

  if (this->d_perturb_map.count(name)) {
    it = this->d_perturb_map.find(name);
    delete std::get<0>(this->d_perturb_map[name]);
    this->d_perturb_map.erase(name);
  } else
    DEBUG_TRACE("This perturb buffer do not exist");

  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int SutraController<Tcomp, Tout>::enable_perturb_voltage(string name) {
  std::lock_guard<std::mutex> lock(this->comp_voltage_mutex);
  current_context->set_active_device(device, 1);
  typename map<string, tuple<CarmaObj<Tcomp> *, int, bool>>::iterator it;

  if (this->d_perturb_map.count(name)) {
    it = this->d_perturb_map.find(name);
    std::get<2>(this->d_perturb_map[name]) = true;
  } else
    DEBUG_TRACE("This perturb buffer do not exist");

  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int SutraController<Tcomp, Tout>::disable_perturb_voltage(string name) {
  std::lock_guard<std::mutex> lock(this->comp_voltage_mutex);
  current_context->set_active_device(device, 1);
  typename map<string, tuple<CarmaObj<Tcomp> *, int, bool>>::iterator it;

  if (this->d_perturb_map.count(name)) {
    it = this->d_perturb_map.find(name);
    std::get<2>(this->d_perturb_map[name]) = false;
  } else
    DEBUG_TRACE("This perturb buffer do not exist");

  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int SutraController<Tcomp, Tout>::reset_perturb_voltage() {
  std::lock_guard<std::mutex> lock(this->comp_voltage_mutex);
  current_context->set_active_device(device, 1);
  typename map<string, tuple<CarmaObj<Tcomp> *, int, bool>>::iterator it;
  it = this->d_perturb_map.begin();
  while (it != this->d_perturb_map.end()) {
    delete std::get<0>(it->second);
    it++;
  }
  this->d_perturb_map.clear();

  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int SutraController<Tcomp, Tout>::command_delay() {
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
int SutraController<Tcomp, Tout>::clip_commands() {
  current_context->set_active_device(device, 1);
  if (this->volt_min < this->volt_max)
    this->d_com_clipped->clip(this->volt_min, this->volt_max);
  else
    DEBUG_TRACE(
        "volt_max value must be greater than volt_min value. Nothing has been done.");
  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int SutraController<Tcomp, Tout>::comp_latency() {
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
int SutraController<Tcomp, Tout>::comp_voltage() {
  std::lock_guard<std::mutex> lock(this->comp_voltage_mutex);
  current_context->set_active_device(device, 1);
  if (this->delay > 0) {
    this->comp_latency();
    this->command_delay();
  } else {
    this->d_com_clipped->copy_from(this->d_com->get_data(),
                                 this->d_com_clipped->get_nb_elements());
  }
  this->add_perturb();
  this->clip_commands();
  convert_to_voltage<Tcomp, Tout>(
      this->d_com_clipped->get_data(), this->d_voltage->get_data(), this->nactu(),
      this->volt_min, this->volt_max, this->val_max,
      this->current_context->get_device(this->device));

  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int SutraController<Tcomp, Tout>::add_perturb() {
  typename map<string, tuple<CarmaObj<Tcomp> *, int, bool>>::iterator it;
  int cpt;
  CarmaObj<Tcomp> *d_perturb;
  int any_perturb = 0;
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
  delete this->streams;

  if (this->d_centroids != nullptr) delete this->d_centroids;
  if (this->d_centroids_padded != nullptr) delete this->d_centroids_padded;
  if (this->d_com_padded != nullptr) delete this->d_com_padded;
  for (int k=0; k < this->d_circular_coms.size() ; k++) {
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
int SutraController<Tcomp, Tout>::set_com(float *com, int nElem) {
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
int SutraController<Tcomp, Tout>::set_delay(float delay) {

  this->delay = delay;
  int floor = (int)delay;
  this->a = Tcomp(delay - floor);
  this->b = Tcomp(1 - a);


  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int SutraController<Tcomp, Tout>::set_volt_min(float volt_min) {
  this->volt_min = volt_min;
  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int SutraController<Tcomp, Tout>::set_volt_max(float volt_max) {
  this->volt_max = volt_max;
  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int SutraController<Tcomp, Tout>::set_val_max(float val_max) {
  this->val_max = Tcomp(val_max);
  return EXIT_SUCCESS;
}

template class SutraController<float, float>;
template class SutraController<float, uint16_t>;

#ifdef CAN_DO_HALF
template class SutraController<half, float>;
template class SutraController<half, uint16_t>;
#endif
// template<class T>
// int SutraController<Tcomp, Tout>::syevd_f(char meth, CarmaObj<float> *d_U,
//                               CarmaHostObj<float> *h_eigenvals) {
//   current_context->set_active_device(device, 1);

//   // Init double arrays
//   const long dims_data[3] = {2, d_U->get_dims()[1], d_U->get_dims()[2]};
//   CarmaObj<double> *d_Udouble =
//       new CarmaObj<double>(current_context, dims_data);
//   const long dims_data2[2] = {1, h_eigenvals->get_dims()[1]};
//   CarmaHostObj<double> *h_eigen_double =
//       new CarmaHostObj<double>(dims_data2, MA_PAGELOCK);

//   // Copy float array in double array
//   float_to_double(d_U->get_data(), d_Udouble->get_data(), d_U->get_nb_elements(),
//                 this->current_context->get_device(device));

//   // Doing syevd<double>
//   carma_magma_syevd<double, 1>(meth, d_Udouble, h_eigen_double);

//   // Reverse copy
//   double_to_float(d_Udouble->get_data(), d_U->get_data(), d_U->get_nb_elements(),
//                 this->current_context->get_device(device));

//   for (int cc = 0; cc < h_eigenvals->get_nb_elements(); cc++) {
//     h_eigenvals->get_data()[cc] = (float)h_eigen_double->get_data()[cc];
//   }

//   delete d_Udouble;
//   delete h_eigen_double;

//   return EXIT_SUCCESS;
// }

// template<class T>
// int SutraController<Tcomp, Tout>::invgen(CarmaObj<float> *d_mat, float
// cond, int job)
// {
//   current_context->set_active_device(device, 1);
//   const long dims_data[3] = {2, d_mat->get_dims()[1], d_mat->get_dims()[2]};
//   CarmaObj<float> *d_U = new CarmaObj<float>(current_context, dims_data);
//   CarmaObj<float> *d_tmp = new CarmaObj<float>(current_context, dims_data);
//   int i;

//   const long dims_data2[2] = {1, d_mat->get_dims()[1]};
//   CarmaObj<float> *d_eigenvals_inv =
//       new CarmaObj<float>(current_context, dims_data2);
//   CarmaHostObj<float> *h_eigenvals =
//       new CarmaHostObj<float>(dims_data2, MA_PAGELOCK);
//   CarmaHostObj<float> *h_eigenvals_inv =
//       new CarmaHostObj<float>(dims_data2, MA_PAGELOCK);

//   d_U->copy(d_mat, 1, 1);
//   carma_magma_syevd<float, 1>('V', d_U, h_eigenvals);

//   // syevd_f('V',d_U,h_eigenvals);
//   if (job == 1) {  // Conditionnement
//     float maxe = h_eigenvals->get_data()[d_mat->get_dims()[1] - 1];

//     for (i = 0; i < d_mat->get_dims()[1]; i++) {
//       if (h_eigenvals->get_data()[i] < maxe / cond)
//         h_eigenvals_inv->get_data()[i] = 0.;
//       else
//         h_eigenvals_inv->get_data()[i] = 1. / h_eigenvals->get_data()[i];
//     }
//   }

//   if (job == 0) {  // Filtre #cond modes
//     for (i = 0; i < d_mat->get_dims()[1]; i++) {
//       if (i < cond)
//         h_eigenvals_inv->get_data()[i] = 0.;
//       else
//         h_eigenvals_inv->get_data()[i] = 1. / h_eigenvals->get_data()[i];
//     }
//   }

//   d_eigenvals_inv->host2device(*h_eigenvals_inv);

//   carma_dgmm(this->cublas_handle(), CUBLAS_SIDE_RIGHT, d_mat->get_dims()[1],
//              d_mat->get_dims()[2], d_U->get_data(), d_mat->get_dims()[1],
//              d_eigenvals_inv->get_data(), 1, d_tmp->get_data(),
//              d_mat->get_dims()[1]);
//   carma_gemm<float>(this->cublas_handle(), 'n', 't', d_mat->get_dims()[1],
//                     d_mat->get_dims()[1], d_mat->get_dims()[2], 1.0f,
//                     d_tmp->get_data(), d_mat->get_dims()[1], d_U->get_data(),
//                     d_mat->get_dims()[1], 0.0f, d_mat->get_data(),
//                     d_mat->get_dims()[1]);

//   delete d_U;
//   delete d_tmp;
//   delete d_eigenvals_inv;
//   delete h_eigenvals;
//   delete h_eigenvals_inv;

//   return EXIT_SUCCESS;
// }
