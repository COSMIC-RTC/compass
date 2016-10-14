#include <sutra_controller.h>
#include <string>

sutra_controller::sutra_controller(carma_context *context,
                                   int            nslope,
                                   int            nactu,
                                   float          delay,
                                   sutra_dms     *dms,
                                   char         **type,
                                   float         *alt,
                                   int            ndm) {
  this->current_context = context;
  this->device          = context->get_activeDevice();

  // current_context->set_activeDevice(device,1);
  this->d_com1 = NULL;
  this->d_com2 = NULL;

  int nstreams = 1; // nvalid/10;

  while (nactu % nstreams != 0) nstreams--;

  std::cerr << "controller uses " << nstreams << " streams" << std::endl;
  streams = new carma_streams(nstreams);

  this->open_loop = 0;
  this->d_perturb = NULL;
  this->cpt_pertu = 0;

  if (delay < 2) this->delay = delay;
  else this->delay = 2.0f;

  int floor = (int)delay;

  if ((floor > 0) && (floor < 2)) {
    this->a = 0;
    this->c = delay - floor;
    this->b = 1 - this->c;
  } else if (floor == 0) {
    this->b = delay;
    this->a = 1 - this->b;
    this->c = 0;
  } else { // Maximum delay is 2
    this->a = 0;
    this->c = 1;
    this->b = 0;
  }

  // DEBUG_TRACE("delay = %f a = %f b = %f c = %f floor =
  // %f",this->delay,this->a,this->b,this->c,floor);

  long dims_data1[2] = { 1, 0 };

  dims_data1[1]     = nslope / 2;
  this->d_subsum    = new carma_obj<float>(context, dims_data1);
  dims_data1[1]     = nslope;
  this->d_centroids = new carma_obj<float>(context, dims_data1);
  dims_data1[1]     = nactu;
  this->d_com       = new carma_obj<float>(context, dims_data1);
  this->d_com1      = new carma_obj<float>(context, dims_data1);

  if (this->delay > 1) {
    this->d_com2 = new carma_obj<float>(context, dims_data1);
  }

  this->d_voltage = new carma_obj<float>(context, dims_data1);

  for (int i = 0; i < ndm; i++) {
    this->d_dmseen.push_back(dms->d_dms[dms->get_inddm(type[i], alt[i])]);
  }
}

int sutra_controller::set_openloop(int open_loop_status) {
  current_context->set_activeDevice(device, 1);
  this->open_loop = open_loop_status;

  if (this->open_loop) {
    carmaSafeCall(
      cudaMemset(this->d_com->getData(), 0.0f,
                 this->nactu() * sizeof(float)));
    carmaSafeCall(
      cudaMemset(this->d_com1->getData(), 0.0f,
                 this->nactu() * sizeof(float)));

    if (this->delay > 1) {
      carmaSafeCall(
        cudaMemset(this->d_com2->getData(), 0.0f,
                   this->nactu() * sizeof(float)));
    }
  }
  return EXIT_SUCCESS;
}

int sutra_controller::set_perturbcom(float *perturb, int N) {
  carma_obj<float> *tmp = nullptr;
  if (this->d_perturb != NULL) tmp = this->d_perturb;
  long dims_data2[3] = { 2, this->nactu(), N };

  current_context->set_activeDevice(device, 1);
  this->d_perturb = new carma_obj<float>(current_context, dims_data2);
  this->d_perturb->host2device(perturb);

  this->cpt_pertu = 0;

  if (tmp != nullptr) delete tmp;

  return EXIT_SUCCESS;
}

int sutra_controller::command_delay() {
  current_context->set_activeDevice(device, 1);

  if (delay > 1) {
    this->d_com2->copy(this->d_com1, 1, 1);
    this->d_com1->copy(this->d_com, 1, 1);
  } else if (delay > 0) this->d_com1->copy(this->d_com, 1, 1);

  return EXIT_SUCCESS;
}

void sutra_controller::clip_voltage(float min, float max) {
  current_context->set_activeDevice(device, 1);
  this->d_com->clip(min, max);
}

int sutra_controller::comp_voltage() {
  current_context->set_activeDevice(device, 1);

  if (this->open_loop)
    carmaSafeCall(
      cudaMemset(this->d_voltage->getData(), 0.0f,
                 this->nactu() * sizeof(float)));
  else if (this->delay > 0) {
    // DEBUG_TRACE("a = %f   b = %f   c = %f",this->a,this->b,this->c);
    carmaSafeCall(
      cudaMemset(this->d_voltage->getData(), 0.0f,
                 this->nactu() * sizeof(float)));
    carma_axpy(cublas_handle(), this->nactu(), this->a, this->d_com->getData(),
               1, this->d_voltage->getData(), 1);
    carma_axpy(cublas_handle(), this->nactu(), this->b, this->d_com1->getData(),
               1, this->d_voltage->getData(), 1);

    if (delay > 1)
      carma_axpy(cublas_handle(), this->nactu(), this->c,
                 this->d_com2->getData(), 1, this->d_voltage->getData(), 1);
  } else this->d_com->copyInto(this->d_voltage->getData(), this->nactu());

  if (this->d_perturb != NULL) { // Apply volt perturbations (circular buffer)
    carma_axpy(cublas_handle(), this->nactu(), 1.0f,
               this->d_perturb->getData(this->cpt_pertu * this->nactu()), 1,
               this->d_voltage->getData(), 1);

    if (this->cpt_pertu < this->d_perturb->getDims()[2] - 1) this->cpt_pertu += 1;
    else this->cpt_pertu = 0;
  }

  return EXIT_SUCCESS;
}

sutra_controller::~sutra_controller() {
  delete this->streams;

  delete this->d_subsum;
  delete this->d_centroids;
  delete this->d_com;
  delete this->d_com1;
  delete this->d_voltage;

  if (this->d_perturb != NULL) delete this->d_perturb;
}

int sutra_controller::syevd_f(char meth, carma_obj<float> *d_U,
                              carma_host_obj<float> *h_eigenvals) {
  current_context->set_activeDevice(device, 1);

  // Init double arrays
  const long dims_data[3]      = { 2, d_U->getDims()[1], d_U->getDims()[2] };
  carma_obj<double> *d_Udouble = new carma_obj<double>(current_context,
                                                       dims_data);
  const long dims_data2[2]               = { 1, h_eigenvals->getDims()[1] };
  carma_host_obj<double> *h_eigen_double = new carma_host_obj<double>(
    dims_data2, MA_PAGELOCK);

  // Copy float array in double array
  floattodouble(d_U->getData(), d_Udouble->getData(), d_U->getNbElem(),
                this->current_context->get_device(device));

  // Doing syevd<double>
  carma_syevd<double, 1>(meth, d_Udouble, h_eigen_double);

  // Reverse copy
  doubletofloat(d_Udouble->getData(), d_U->getData(), d_U->getNbElem(),
                this->current_context->get_device(device));

  for (int cc = 0; cc < h_eigenvals->getNbElem(); cc++) {
    h_eigenvals->getData()[cc] = (float)h_eigen_double->getData()[cc];
  }

  delete d_Udouble;
  delete h_eigen_double;

  return EXIT_SUCCESS;
}

int sutra_controller::invgen(carma_obj<float> *d_mat, float cond, int job) {
  current_context->set_activeDevice(device, 1);
  const long dims_data[3] = { 2, d_mat->getDims()[1], d_mat->getDims()[2] };
  carma_obj<float> *d_U   = new carma_obj<float>(current_context, dims_data);
  carma_obj<float> *d_tmp = new carma_obj<float>(current_context, dims_data);
  int i;

  const long dims_data2[2]          = { 1, d_mat->getDims()[1] };
  carma_obj<float> *d_eigenvals_inv = new carma_obj<float>(current_context,
                                                           dims_data2);
  carma_host_obj<float> *h_eigenvals = new carma_host_obj<float>(dims_data2,
                                                                 MA_PAGELOCK);
  carma_host_obj<float> *h_eigenvals_inv = new carma_host_obj<float>(
    dims_data2, MA_PAGELOCK);

  d_U->copy(d_mat, 1, 1);
  carma_syevd<float, 1>('V', d_U, h_eigenvals);

  // syevd_f('V',d_U,h_eigenvals);
  if (job == 1) { // Conditionnement
    float maxe = h_eigenvals->getData()[d_mat->getDims()[1] - 1];

    for (i = 0; i < d_mat->getDims()[1]; i++) {
      if (h_eigenvals->getData()[i] < maxe / cond) h_eigenvals_inv->getData()[i] = 0.;
      else h_eigenvals_inv->getData()[i] = 1. / h_eigenvals->getData()[i];
    }
  }

  if (job == 0) { // Filtre #cond modes
    for (i = 0; i < d_mat->getDims()[1]; i++) {
      if (i < cond) h_eigenvals_inv->getData()[i] = 0.;
      else h_eigenvals_inv->getData()[i] = 1. / h_eigenvals->getData()[i];
    }
  }

  d_eigenvals_inv->host2device(*h_eigenvals_inv);

  carma_dgmm(cublas_handle(), CUBLAS_SIDE_RIGHT, d_mat->getDims()[1],
             d_mat->getDims()[2], d_U->getData(), d_mat->getDims()[1],
             d_eigenvals_inv->getData(), 1, d_tmp->getData(),
             d_mat->getDims()[1]);
  carma_gemm<float>(cublas_handle(), 'n', 't', d_mat->getDims()[1],
                    d_mat->getDims()[1], d_mat->getDims()[2], 1.0f,
                    d_tmp->getData(), d_mat->getDims()[1], d_U->getData(),
                    d_mat->getDims()[1], 0.0f, d_mat->getData(),
                    d_mat->getDims()[1]);

  delete d_U;
  delete d_tmp;
  delete d_eigenvals_inv;
  delete h_eigenvals;
  delete h_eigenvals_inv;

  return EXIT_SUCCESS;
}
