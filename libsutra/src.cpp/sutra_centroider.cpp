#include <sutra_centroider.h>

sutra_centroider::sutra_centroider(carma_context *context, sutra_wfs *wfs,
                                   long nvalid, float offset, float scale,
                                   int device) {
  this->current_context = context;
  this->device = device;
  context->set_activeDevice(device, 1);
  this->wfs = wfs;

  this->nvalid = nvalid;
  this->offset = offset;
  this->scale = scale;
  this->nslopes = 0;

  this->d_bincube = nullptr;
  this->d_subsum = nullptr;
  this->d_img = nullptr;
  this->d_validx = nullptr;
  this->d_validy = nullptr;
  this->d_dark = nullptr;
  this->d_flat = nullptr;
}

sutra_centroider::~sutra_centroider() {
  if (this->d_bincube != nullptr) delete this->d_bincube;
  if (this->d_img != nullptr) delete this->d_img;
  if (this->d_validx != nullptr) delete this->d_validx;
  if (this->d_validy != nullptr) delete this->d_validy;
  if (this->d_dark != nullptr) delete this->d_dark;
  if (this->d_flat != nullptr) delete this->d_flat;
}

int sutra_centroider::set_scale(float scale) {
  this->scale = scale;
  return EXIT_SUCCESS;
}

int sutra_centroider::set_dark(float *dark, int n) {
  current_context->set_activeDevice(device, 1);
  if (this->d_dark == nullptr) {
    long dims_data2[3] = {2, n, n};
    this->d_dark = new carma_obj<float>(current_context, dims_data2);
  }
  this->d_dark->host2device(dark);
  return EXIT_SUCCESS;
}

int sutra_centroider::set_flat(float *flat, int n) {
  current_context->set_activeDevice(device, 1);
  if (this->d_flat == nullptr) {
    long dims_data2[3] = {2, n, n};
    this->d_flat = new carma_obj<float>(current_context, dims_data2);
  }
  this->d_flat->host2device(flat);
  return EXIT_SUCCESS;
}

int sutra_centroider::calibrate_img() {
  carma_obj<float> *img;

  if (this->d_img != nullptr)
    img = this->d_img;
  else
    img = this->d_bincube;

  if (this->d_dark != nullptr) img->axpy(-1.f, this->d_dark, 1, 1);

  if (this->d_flat != nullptr)
    mult_vect(img->getData(), this->d_flat->getData(), img->getNbElem(),
              current_context->get_device(device));

  return EXIT_SUCCESS;
}
int sutra_centroider::load_img(float *img, int n) {
  current_context->set_activeDevice(device, 1);
  if (this->d_img == nullptr) {
    long dims_data2[3] = {2, n, n};
    this->d_img = new carma_obj<float>(current_context, dims_data2);
  }
  this->d_img->host2device(img);
  return EXIT_SUCCESS;
}

int sutra_centroider::load_img_gpu(float *img, int n) {
  current_context->set_activeDevice(device, 1);
  if (this->d_img == nullptr) {
    long dims_data2[3] = {2, n, n};
    this->d_img = new carma_obj<float>(current_context, dims_data2);
  }
  this->d_img->copyFrom(img, n * n);
  return EXIT_SUCCESS;
}

int sutra_centroider::load_pyrimg(float *img, int n) {
  current_context->set_activeDevice(device, 1);
  if (this->d_bincube == nullptr) {
    long dims_data2[3] = {2, n, n};
    this->d_bincube = new carma_obj<float>(current_context, dims_data2);
  }
  this->d_bincube->host2device(img);
  return EXIT_SUCCESS;
}

int sutra_centroider::fill_bincube(int npix) {
  current_context->set_activeDevice(device, 1);
  if (this->d_bincube == nullptr) {
    long dims_data3[4] = {3, npix, npix, this->nvalid};
    this->d_bincube = new carma_obj<float>(current_context, dims_data3);
  }
  int nxsub = this->d_img->getDims(1) / npix;
  fillbincube(this->d_img->getData(), this->d_bincube->getData(), npix,
              this->nvalid, nxsub, this->d_validx->getData(),
              this->d_validy->getData(),
              this->current_context->get_device(this->device));

  return EXIT_SUCCESS;
}

int sutra_centroider::load_validpos(int *ivalid, int *jvalid, int N) {
  current_context->set_activeDevice(device, 1);
  if (this->d_validx == nullptr) {
    long dims_data[2] = {1, N};
    this->d_validx = new carma_obj<int>(current_context, dims_data);
    this->d_validy = new carma_obj<int>(current_context, dims_data);
  }

  this->d_validx->host2device(ivalid);
  this->d_validy->host2device(jvalid);

  return EXIT_SUCCESS;
}
