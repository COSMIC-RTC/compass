#include <sutra_controller.h>
#include <string>

template <class T>
sutra_controller<T>::sutra_controller(carma_context *context, int nvalid,
                                      int nslope, int nactu, T delay,
                                      sutra_dms *dms, int *idx_dms, int ndm) {
  this->current_context = context;
  this->device = context->get_activeDevice();

  // current_context->set_activeDevice(device,1);
  this->d_com1 = NULL;
  this->d_com2 = NULL;

  int nstreams = 1;  // nvalid/10;

  while (nactu % nstreams != 0) nstreams--;

  std::cerr << "controller uses " << nstreams << " streams" << std::endl;
  streams = new carma_streams(nstreams);

  this->open_loop = 0;

  if (delay < 2)
    this->delay = delay;
  else
    this->delay = 2.0f;

  int floor = (int)delay;

  if ((floor > 0) && (floor < 2)) {
    this->a = 0.f;
    this->c = delay - floor;
    this->b = 1 - this->c;
  } else if (floor == 0) {
    this->b = delay;
    this->a = 1 - this->b;
    this->c = 0.f;
  } else {  // Maximum delay is 2
    this->a = 0.f;
    this->c = 1.f;
    this->b = 0.f;
  }

  // DEBUG_TRACE("delay = %f a = %f b = %f c = %f floor =
  // %f",this->delay,this->a,this->b,this->c,floor);

  long dims_data1[2] = {1, 0};

  dims_data1[1] = nslope;
  this->d_centroids = new carma_obj<T>(context, dims_data1);

  dims_data1[1] = nactu;
  this->d_com = new carma_obj<T>(context, dims_data1);
  this->d_com1 = new carma_obj<T>(context, dims_data1);

  if (this->delay > 1) {
    this->d_com2 = new carma_obj<T>(context, dims_data1);
  }

  this->d_voltage = new carma_obj<T>(context, dims_data1);

  for (int i = 0; i < ndm; i++) {
    this->d_dmseen.push_back(dms->d_dms[idx_dms[i]]);
  }
}

template <class T>
int sutra_controller<T>::set_openloop(int open_loop_status, bool rst) {
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

template <class T>
int sutra_controller<T>::add_perturb_voltage(string name, T *perturb, int N) {
  std::lock_guard<std::mutex> lock(this->comp_voltage_mutex);
  current_context->set_activeDevice(device, 1);

  if (this->d_perturb_map.count(name) < 1) {
    long dims_data2[3] = {2, N, this->nactu()};
    int cpt = 0;
    this->d_perturb_map[name] = std::make_tuple(
        new carma_obj<T>(current_context, dims_data2, perturb), cpt, true);
  } else
    DEBUG_TRACE("This perturb buffer already exists");

  return EXIT_SUCCESS;
}

template <class T>
int sutra_controller<T>::remove_perturb_voltage(string name) {
  std::lock_guard<std::mutex> lock(this->comp_voltage_mutex);
  current_context->set_activeDevice(device, 1);
  typename map<string, tuple<carma_obj<T> *, int, bool>>::iterator it;

  if (this->d_perturb_map.count(name)) {
    it = this->d_perturb_map.find(name);
    delete std::get<0>(this->d_perturb_map[name]);
    this->d_perturb_map.erase(name);
  } else
    DEBUG_TRACE("This perturb buffer do not exist");

  return EXIT_SUCCESS;
}

template <class T>
int sutra_controller<T>::enable_perturb_voltage(string name) {
  std::lock_guard<std::mutex> lock(this->comp_voltage_mutex);
  current_context->set_activeDevice(device, 1);
  typename map<string, tuple<carma_obj<T> *, int, bool>>::iterator it;

  if (this->d_perturb_map.count(name)) {
    it = this->d_perturb_map.find(name);
    std::get<2>(this->d_perturb_map[name]) = true;
  } else
    DEBUG_TRACE("This perturb buffer do not exist");

  return EXIT_SUCCESS;
}

template <class T>
int sutra_controller<T>::disable_perturb_voltage(string name) {
  std::lock_guard<std::mutex> lock(this->comp_voltage_mutex);
  current_context->set_activeDevice(device, 1);
  typename map<string, tuple<carma_obj<T> *, int, bool>>::iterator it;

  if (this->d_perturb_map.count(name)) {
    it = this->d_perturb_map.find(name);
    std::get<2>(this->d_perturb_map[name]) = false;
  } else
    DEBUG_TRACE("This perturb buffer do not exist");

  return EXIT_SUCCESS;
}

template <class T>
int sutra_controller<T>::reset_perturb_voltage() {
  std::lock_guard<std::mutex> lock(this->comp_voltage_mutex);
  current_context->set_activeDevice(device, 1);
  typename map<string, tuple<carma_obj<T> *, int, bool>>::iterator it;
  it = this->d_perturb_map.begin();
  while (it != this->d_perturb_map.end()) {
    delete std::get<0>(it->second);
    it++;
  }
  this->d_perturb_map.clear();

  return EXIT_SUCCESS;
}

template <class T>
int sutra_controller<T>::command_delay() {
  current_context->set_activeDevice(device, 1);

  if (delay > 1) {
    this->d_com2->copy(this->d_com1, 1, 1);
    this->d_com1->copy(this->d_com, 1, 1);
  } else if (delay > 0)
    this->d_com1->copy(this->d_com, 1, 1);

  return EXIT_SUCCESS;
}

template <class T>
void sutra_controller<T>::clip_voltage(T min, T max) {
  current_context->set_activeDevice(device, 1);
  if (min < max)
    this->d_com->clip(min, max);
  else
    DEBUG_TRACE(
        "max value must be greater than min value. Nothing has been done.");
}

template <class T>
int sutra_controller<T>::comp_voltage() {
  std::lock_guard<std::mutex> lock(this->comp_voltage_mutex);
  current_context->set_activeDevice(device, 1);

  this->d_voltage->reset();

  if (!this->open_loop) {
    if (this->delay > 0) {
      this->d_voltage->axpy(this->a, this->d_com, 1, 1);
      this->d_voltage->axpy(this->b, this->d_com1, 1, 1);

      if (delay > 1) this->d_voltage->axpy(this->b, this->d_com2, 1, 1);

    } else {
      this->d_com->copyInto(this->d_voltage->getData(), this->nactu());
    }
  } else {
    this->d_com->copyInto(this->d_voltage->getData(), this->nactu());
  }

  add_perturb();
  return EXIT_SUCCESS;
}

template <class T>
int sutra_controller<T>::add_perturb() {
  typename map<string, tuple<carma_obj<T> *, int, bool>>::iterator it;
  int cpt;
  carma_obj<T> *d_perturb;
  for (it = this->d_perturb_map.begin(); it != this->d_perturb_map.end();
       ++it) {
    if (std::get<2>(it->second)) {
      cpt = std::get<1>(it->second);
      d_perturb = std::get<0>(it->second);
      this->d_voltage->axpy(1.0f, d_perturb, 1, 1, cpt);
      if (cpt < d_perturb->getDims(1) - 1)
        std::get<1>(it->second) = cpt + 1;
      else
        std::get<1>(it->second) = 0;
    }
  }

  return EXIT_SUCCESS;
}

template <class T>
sutra_controller<T>::~sutra_controller() {
  delete this->streams;

  delete this->d_centroids;
  delete this->d_com;
  delete this->d_com1;
  delete this->d_voltage;
  this->reset_perturb_voltage();
}

template <class T>
int sutra_controller<T>::set_com(T *com, int nElem) {
  current_context->set_activeDevice(device, 1);
  if (nElem == this->d_com->getNbElem())
    this->d_com->host2device(com);
  else
    DEBUG_TRACE("Wrong dimension of com");

  return EXIT_SUCCESS;
}

template class sutra_controller<float>;

#ifdef CAN_DO_HALF
template class sutra_controller<half>;
#endif
// template<class T>
// int sutra_controller<T>::syevd_f(char meth, carma_obj<float> *d_U,
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
// int sutra_controller<T>::invgen(carma_obj<float> *d_mat, float cond, int job)
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
