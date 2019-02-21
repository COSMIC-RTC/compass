#include <sutra_controller.h>
#include <string>

template <typename Tcomp, typename Tout>
sutra_controller<Tcomp, Tout>::sutra_controller(carma_context *context,
                                                int nvalid, int nslope,
                                                int nactu, float delay,
                                                sutra_dms *dms, int *idx_dms,
                                                int ndm) {
  this->current_context = context;
  this->device = context->get_activeDevice();

  this->nactus = nactu;
  this->nslopes = nslope;
  // current_context->set_activeDevice(device,1);
  this->d_com1 = nullptr;
  this->d_com2 = nullptr;
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
  if (delay < 2)
    this->delay = Tcomp(delay);
  else
    this->delay = Tcomp(2.0f);

  int floor = (int)delay;

  if ((floor > 0) && (floor < 2)) {
    this->a = Tcomp(0.f);
    this->c = Tcomp(delay - floor);
    this->b = Tcomp(1.f) - this->c;
  } else if (floor == 0) {
    this->b = this->delay;
    this->a = Tcomp(1.f) - this->b;
    this->c = Tcomp(0.f);
  } else {  // Maximum delay is 2
    this->a = Tcomp(0.f);
    this->c = Tcomp(1.f);
    this->b = Tcomp(0.f);
  }

  // DEBUG_TRACE("delay = %f a = %f b = %f c = %f floor =
  // %f",this->delay,this->a,this->b,this->c,floor);

  long dims_data1[2] = {1, 0};

  dims_data1[1] = nslope;
  this->d_centroids = new carma_obj<Tcomp>(context, dims_data1);

  dims_data1[1] = nactu;
  this->d_com = new carma_obj<Tcomp>(context, dims_data1);
  this->d_com1 = new carma_obj<Tcomp>(context, dims_data1);
  this->d_comClipped = new carma_obj<Tcomp>(context, dims_data1);

  if (this->delay > 1) {
    this->d_com2 = new carma_obj<Tcomp>(context, dims_data1);
  }

  this->init_voltage();

  for (int i = 0; i < ndm; i++) {
    this->d_dmseen.push_back(dms->d_dms[idx_dms[i]]);
  }

  if (std::is_same<Tcomp, half>::value) {
    while (nslope % 8 != 0)
      nslope++;  // Watching for multiple of 8 to get max perf with tensor cores
    while (nactu % 8 != 0) nactu++;
    dims_data1[1] = nslope;
    this->d_centroidsPadded = new carma_obj<Tcomp>(context, dims_data1);
    this->d_centroids->swapPtr(this->d_centroidsPadded->getData());
    dims_data1[1] = nactu;
    this->d_comPadded = new carma_obj<Tcomp>(context, dims_data1);
    this->d_com->swapPtr(this->d_comPadded->getData());
  }
}

template <typename Tcomp, typename Tout>
int sutra_controller<Tcomp, Tout>::set_openloop(int open_loop_status,
                                                bool rst) {
  current_context->set_activeDevice(device, 1);
  this->open_loop = open_loop_status;

  if (this->open_loop && rst) {
    this->d_com->reset();
    this->d_com1->reset();

    if (this->delay > 1) {
      this->d_com2->reset();
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

  if (delay > 1) {
    this->d_com2->copyFrom(this->d_com1->getData(), this->d_com2->getNbElem());
    this->d_com1->copyFrom(this->d_com->getData(), this->d_com1->getNbElem());
  } else if (delay > 0)
    this->d_com1->copyFrom(this->d_com->getData(), this->d_com1->getNbElem());

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
    this->d_comClipped->axpy(this->a, this->d_com, 1, 1);
    this->d_comClipped->axpy(this->b, this->d_com1, 1, 1);

    if (delay > 1) this->d_comClipped->axpy(this->c, this->d_com2, 1, 1);
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
                                 this->d_com->getNbElem());
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
  if (this->d_com != nullptr) delete this->d_com;
  if (this->d_com1 != nullptr) delete this->d_com1;
  if (this->d_com2 != nullptr) delete this->d_com2;
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
                                 this->d_com->getNbElem());
  } else
    DEBUG_TRACE("Wrong dimension of com");

  return EXIT_SUCCESS;
}

template <typename Tcomp, typename Tout>
int sutra_controller<Tcomp, Tout>::set_delay(float delay) {
  this->delay = Tcomp(delay);
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
//   carma_syevd<double, 1>(meth, d_Udouble, h_eigen_double);

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
//   carma_syevd<float, 1>('V', d_U, h_eigenvals);

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
