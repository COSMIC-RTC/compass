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
//! \class     sutra_controller
//! \brief     this class provides the controller features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.3.2
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#include <sutra_controller.h>
#include <string>

template <typename Tcomp, typename Tout>
sutra_controller<Tcomp, Tout>::sutra_controller(carma_context *context,
                                                int nvalid, int nslope,
                                                int nactu, float delay,
                                                sutra_dms *dms, int *idx_dms,
                                                int ndm, int *idx_centro, int ncentro) {
  this->current_context = context;
  this->device = context->get_activeDevice();

  this->nactus = nactu;
  this->nslopes = nslope;
  // current_context->set_activeDevice(device,1);
  this->d_com1 = nullptr;
  this->d_comPadded = nullptr;
  this->d_centroidsPadded = nullptr;

  int nstreams = 1;  // nvalid/10;

  while (nactu % nstreams != 0) nstreams--;

  std::cerr << "controller uses " << nstreams << " streams" << std::endl;
  streams = new carma_streams(nstreams);

  this->open_loop = 0;
  if (std::is_same<Tout, uint16_t>::value) {
    this->Vmin = -1.0f;
    this->Vmax = 1.0f;
  } else {
    this->Vmin = -1.e6;
    this->Vmax = 1.e6;
  }
  this->valMax = Tout(65535);

  this->set_delay(delay);

  long dims_data1[2] = {1, 0};

  dims_data1[1] = nslope;
  this->d_centroids = new carma_obj<Tcomp>(context, dims_data1);

  dims_data1[1] = nactu;
  this->d_com = new carma_obj<Tcomp>(context, dims_data1);
  this->d_circularComs.push_front(d_com);
  this->d_com1 = new carma_obj<Tcomp>(context, dims_data1);
  if (this->delay > 0) {
    this->d_circularComs.push_front(d_com1);
    while (this->d_circularComs.size() <= int(this->delay) + 1) {
      this->d_circularComs.push_front(new carma_obj<Tcomp>(context, dims_data1));
    }
  } 
  this->d_comClipped = new carma_obj<Tcomp>(context, dims_data1);

  this->init_voltage();

  for (int i = 0; i < ndm; i++) {
    this->d_dmseen.push_back(dms->d_dms[idx_dms[i]]);
  }

  for (int i = 0; i < ncentro; i++) {
    this->centro_idx.push_back(idx_centro[i]);
  }

}

template <typename Tcomp, typename Tout>
int sutra_controller<Tcomp, Tout>::set_gain(float gain) {
  this->gain = Tcomp(gain);
  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int sutra_controller<Tcomp, Tout>::set_openloop(int open_loop_status,
                                                bool rst) {
  current_context->set_activeDevice(device, 1);
  this->open_loop = open_loop_status;

  if (this->open_loop && rst) {
    for (auto cobj : this->d_circularComs) {
      cobj->reset();
    }
  }
  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int sutra_controller<Tcomp, Tout>::add_perturb_voltage(string name,
                                                       float *perturb, int N) {
  std::lock_guard<std::mutex> lock(this->comp_voltage_mutex);
  current_context->set_activeDevice(device, 1);

  if (this->d_perturb_map.count(name) < 1) {
    long dims_data2[3] = {2, N, this->nactu()};
    int cpt = 0;
    carma_obj<Tcomp> *d_perturb =
        new carma_obj<Tcomp>(current_context, dims_data2);
    d_perturb->host2device(perturb);
    this->d_perturb_map[name] = std::make_tuple(d_perturb, cpt, true);
  } else
    DEBUG_TRACE("This perturb buffer already exists");

  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int sutra_controller<Tcomp, Tout>::set_perturb_voltage(string name,
                                                       float *perturb, int N) {
  current_context->set_activeDevice(device, 1);

  if (this->d_perturb_map.count(name)) {
    remove_perturb_voltage(name);
  }
  add_perturb_voltage(name, perturb, N);

  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int sutra_controller<Tcomp, Tout>::remove_perturb_voltage(string name) {
  std::lock_guard<std::mutex> lock(this->comp_voltage_mutex);
  current_context->set_activeDevice(device, 1);
  typename map<string, tuple<carma_obj<Tcomp> *, int, bool>>::iterator it;

  if (this->d_perturb_map.count(name)) {
    it = this->d_perturb_map.find(name);
    delete std::get<0>(this->d_perturb_map[name]);
    this->d_perturb_map.erase(name);
  } else
    DEBUG_TRACE("This perturb buffer do not exist");

  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int sutra_controller<Tcomp, Tout>::enable_perturb_voltage(string name) {
  std::lock_guard<std::mutex> lock(this->comp_voltage_mutex);
  current_context->set_activeDevice(device, 1);
  typename map<string, tuple<carma_obj<Tcomp> *, int, bool>>::iterator it;

  if (this->d_perturb_map.count(name)) {
    it = this->d_perturb_map.find(name);
    std::get<2>(this->d_perturb_map[name]) = true;
  } else
    DEBUG_TRACE("This perturb buffer do not exist");

  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int sutra_controller<Tcomp, Tout>::disable_perturb_voltage(string name) {
  std::lock_guard<std::mutex> lock(this->comp_voltage_mutex);
  current_context->set_activeDevice(device, 1);
  typename map<string, tuple<carma_obj<Tcomp> *, int, bool>>::iterator it;

  if (this->d_perturb_map.count(name)) {
    it = this->d_perturb_map.find(name);
    std::get<2>(this->d_perturb_map[name]) = false;
  } else
    DEBUG_TRACE("This perturb buffer do not exist");

  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int sutra_controller<Tcomp, Tout>::reset_perturb_voltage() {
  std::lock_guard<std::mutex> lock(this->comp_voltage_mutex);
  current_context->set_activeDevice(device, 1);
  typename map<string, tuple<carma_obj<Tcomp> *, int, bool>>::iterator it;
  it = this->d_perturb_map.begin();
  while (it != this->d_perturb_map.end()) {
    delete std::get<0>(it->second);
    it++;
  }
  this->d_perturb_map.clear();

  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int sutra_controller<Tcomp, Tout>::command_delay() {
  current_context->set_activeDevice(device, 1);

  if (delay > 0) {
    carma_obj<Tcomp> *tmp = std::move(this->d_circularComs.front());
    tmp->copy(this->d_com, 1, 1);
    this->d_circularComs.pop_front();
    this->d_circularComs.push_back(tmp);
    this->d_com = this->d_circularComs.back();
    this->d_com1 = this->d_circularComs[this->d_circularComs.size() - 2];
  }

  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int sutra_controller<Tcomp, Tout>::clip_commands() {
  current_context->set_activeDevice(device, 1);
  if (this->Vmin < this->Vmax)
    this->d_comClipped->clip(this->Vmin, this->Vmax);
  else
    DEBUG_TRACE(
        "Vmax value must be greater than Vmin value. Nothing has been done.");
  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int sutra_controller<Tcomp, Tout>::comp_latency() {
  // Command increment = a * d[k] + b*d[k-1] + c*d[k-2] + perturb
  this->current_context->set_activeDevice(this->device, 1);
  if (this->delay > 0) {
    this->d_comClipped->reset();
    this->d_comClipped->axpy(this->a, this->d_circularComs[0], 1, 1);
    this->d_comClipped->axpy(this->b, this->d_circularComs[1], 1, 1);
  }

  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int sutra_controller<Tcomp, Tout>::comp_voltage() {
  std::lock_guard<std::mutex> lock(this->comp_voltage_mutex);
  current_context->set_activeDevice(device, 1);
  if (this->delay > 0) {
    this->comp_latency();
    this->command_delay();
  } else {
    this->d_comClipped->copyFrom(this->d_com->getData(),
                                 this->d_comClipped->getNbElem());
  }
  this->add_perturb();
  this->clip_commands();
  convertToVoltage<Tcomp, Tout>(
      this->d_comClipped->getData(), this->d_voltage->getData(), this->nactu(),
      this->Vmin, this->Vmax, this->valMax,
      this->current_context->get_device(this->device));

  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int sutra_controller<Tcomp, Tout>::add_perturb() {
  typename map<string, tuple<carma_obj<Tcomp> *, int, bool>>::iterator it;
  int cpt;
  carma_obj<Tcomp> *d_perturb;
  int any_perturb = 0;
  for (it = this->d_perturb_map.begin(); it != this->d_perturb_map.end();
       ++it) {
    if (std::get<2>(it->second)) {
      cpt = std::get<1>(it->second);
      d_perturb = std::get<0>(it->second);
      this->d_comClipped->axpy(1.0f, d_perturb, d_perturb->getDims(1), 1, cpt);
      if (cpt < d_perturb->getDims(1) - 1)
        std::get<1>(it->second) = cpt + 1;
      else
        std::get<1>(it->second) = 0;
    }
  }

  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
sutra_controller<Tcomp, Tout>::~sutra_controller() {
  delete this->streams;

  if (this->d_centroids != nullptr) delete this->d_centroids;
  if (this->d_centroidsPadded != nullptr) delete this->d_centroidsPadded;
  if (this->d_comPadded != nullptr) delete this->d_comPadded;
  for (int k=0; k < this->d_circularComs.size() ; k++) {
    delete this->d_circularComs[k];
  }
  this->d_circularComs.clear();
  if (this->d_comClipped != nullptr &&
      (void *)this->d_comClipped != (void *)this->d_com)
    delete this->d_comClipped;
  if (this->d_voltage != nullptr &&
      (void *)this->d_voltage != (void *)this->d_comClipped)
    delete this->d_voltage;
  this->reset_perturb_voltage();

}

template <typename Tcomp, typename Tout>
int sutra_controller<Tcomp, Tout>::set_com(float *com, int nElem) {
  current_context->set_activeDevice(device, 1);
  if (nElem == this->d_com->getNbElem()) {
    this->d_com->host2device(com);
    this->d_comClipped->copyFrom(this->d_com->getData(),
                                 this->d_comClipped->getNbElem());
  } else
    DEBUG_TRACE("Wrong dimension of com");

  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int sutra_controller<Tcomp, Tout>::set_delay(float delay) {

  this->delay = delay;
  int floor = (int)delay;
  this->a = Tcomp(delay - floor);
  this->b = Tcomp(1 - a);
  

  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int sutra_controller<Tcomp, Tout>::set_Vmin(float Vmin) {
  this->Vmin = Vmin;
  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int sutra_controller<Tcomp, Tout>::set_Vmax(float Vmax) {
  this->Vmax = Vmax;
  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int sutra_controller<Tcomp, Tout>::set_valMax(float valMax) {
  this->valMax = Tcomp(valMax);
  return EXIT_SUCCESS;
}

template class sutra_controller<float, float>;
template class sutra_controller<float, uint16_t>;

#ifdef CAN_DO_HALF
template class sutra_controller<half, float>;
template class sutra_controller<half, uint16_t>;
#endif
// template<class T>
// int sutra_controller<Tcomp, Tout>::syevd_f(char meth, carma_obj<float> *d_U,
//                               carma_host_obj<float> *h_eigenvals) {
//   current_context->set_activeDevice(device, 1);

//   // Init double arrays
//   const long dims_data[3] = {2, d_U->getDims()[1], d_U->getDims()[2]};
//   carma_obj<double> *d_Udouble =
//       new carma_obj<double>(current_context, dims_data);
//   const long dims_data2[2] = {1, h_eigenvals->getDims()[1]};
//   carma_host_obj<double> *h_eigen_double =
//       new carma_host_obj<double>(dims_data2, MA_PAGELOCK);

//   // Copy float array in double array
//   floattodouble(d_U->getData(), d_Udouble->getData(), d_U->getNbElem(),
//                 this->current_context->get_device(device));

//   // Doing syevd<double>
//   carma_magma_syevd<double, 1>(meth, d_Udouble, h_eigen_double);

//   // Reverse copy
//   doubletofloat(d_Udouble->getData(), d_U->getData(), d_U->getNbElem(),
//                 this->current_context->get_device(device));

//   for (int cc = 0; cc < h_eigenvals->getNbElem(); cc++) {
//     h_eigenvals->getData()[cc] = (float)h_eigen_double->getData()[cc];
//   }

//   delete d_Udouble;
//   delete h_eigen_double;

//   return EXIT_SUCCESS;
// }

// template<class T>
// int sutra_controller<Tcomp, Tout>::invgen(carma_obj<float> *d_mat, float
// cond, int job)
// {
//   current_context->set_activeDevice(device, 1);
//   const long dims_data[3] = {2, d_mat->getDims()[1], d_mat->getDims()[2]};
//   carma_obj<float> *d_U = new carma_obj<float>(current_context, dims_data);
//   carma_obj<float> *d_tmp = new carma_obj<float>(current_context, dims_data);
//   int i;

//   const long dims_data2[2] = {1, d_mat->getDims()[1]};
//   carma_obj<float> *d_eigenvals_inv =
//       new carma_obj<float>(current_context, dims_data2);
//   carma_host_obj<float> *h_eigenvals =
//       new carma_host_obj<float>(dims_data2, MA_PAGELOCK);
//   carma_host_obj<float> *h_eigenvals_inv =
//       new carma_host_obj<float>(dims_data2, MA_PAGELOCK);

//   d_U->copy(d_mat, 1, 1);
//   carma_magma_syevd<float, 1>('V', d_U, h_eigenvals);

//   // syevd_f('V',d_U,h_eigenvals);
//   if (job == 1) {  // Conditionnement
//     float maxe = h_eigenvals->getData()[d_mat->getDims()[1] - 1];

//     for (i = 0; i < d_mat->getDims()[1]; i++) {
//       if (h_eigenvals->getData()[i] < maxe / cond)
//         h_eigenvals_inv->getData()[i] = 0.;
//       else
//         h_eigenvals_inv->getData()[i] = 1. / h_eigenvals->getData()[i];
//     }
//   }

//   if (job == 0) {  // Filtre #cond modes
//     for (i = 0; i < d_mat->getDims()[1]; i++) {
//       if (i < cond)
//         h_eigenvals_inv->getData()[i] = 0.;
//       else
//         h_eigenvals_inv->getData()[i] = 1. / h_eigenvals->getData()[i];
//     }
//   }

//   d_eigenvals_inv->host2device(*h_eigenvals_inv);

//   carma_dgmm(this->cublas_handle(), CUBLAS_SIDE_RIGHT, d_mat->getDims()[1],
//              d_mat->getDims()[2], d_U->getData(), d_mat->getDims()[1],
//              d_eigenvals_inv->getData(), 1, d_tmp->getData(),
//              d_mat->getDims()[1]);
//   carma_gemm<float>(this->cublas_handle(), 'n', 't', d_mat->getDims()[1],
//                     d_mat->getDims()[1], d_mat->getDims()[2], 1.0f,
//                     d_tmp->getData(), d_mat->getDims()[1], d_U->getData(),
//                     d_mat->getDims()[1], 0.0f, d_mat->getData(),
//                     d_mat->getDims()[1]);

//   delete d_U;
//   delete d_tmp;
//   delete d_eigenvals_inv;
//   delete h_eigenvals;
//   delete h_eigenvals_inv;

//   return EXIT_SUCCESS;
// }
