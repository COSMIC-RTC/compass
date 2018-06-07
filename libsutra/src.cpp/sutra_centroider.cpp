#include <sutra_centroider.h>

sutra_centroider::sutra_centroider(carma_context *context,
                                   sutra_sensors *sensors, int nwfs,
                                   long nvalid, float offset, float scale,
                                   int device) {
  this->current_context = context;
  this->device = device;
  context->set_activeDevice(device, 1);
  if (sensors != nullptr)
    this->wfs = sensors->d_wfs[nwfs];
  else
    this->wfs = nullptr;
  this->nwfs = nwfs;
  this->nvalid = nvalid;
  this->offset = offset;
  this->scale = scale;

  this->d_bincube = nullptr;
  this->d_img = nullptr;
  this->d_validx = nullptr;
  this->d_validy = nullptr;
}

sutra_centroider::~sutra_centroider() {
  if (this->d_bincube != nullptr) delete this->d_bincube;
  if (this->d_img != nullptr) delete this->d_img;
  if (this->d_validx != nullptr) delete this->d_validx;
  if (this->d_validy != nullptr) delete this->d_validy;
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
