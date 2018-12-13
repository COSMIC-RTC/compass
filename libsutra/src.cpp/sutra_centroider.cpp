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
  this->npix = 0;

  long dims_data[2] = {1, this->nvalid};
  this->d_intensities = new carma_obj<float>(current_context, dims_data);
  this->d_intensities->reset();

  this->d_centroids_ref = nullptr;
  this->d_img = nullptr;
  this->d_img_raw = nullptr;
  this->d_validx = nullptr;
  this->d_validy = nullptr;
  this->d_dark = nullptr;
  this->d_flat = nullptr;
}

sutra_centroider::~sutra_centroider() {
  if (this->d_intensities != nullptr) delete this->d_intensities;
  if (this->d_centroids_ref != nullptr) delete this->d_centroids_ref;
  if (this->d_img != nullptr) delete this->d_img;
  if (this->d_img_raw != nullptr) delete this->d_img_raw;
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

int sutra_centroider::calibrate_img(bool save_raw) {
  current_context->set_activeDevice(device, 1);

  if (this->d_img == nullptr) {
    std::cout << "Image not initialized\n" << std::endl;
    return EXIT_FAILURE;
  }

  if (save_raw) {
    if (this->d_img_raw != nullptr)
      this->d_img_raw = new carma_obj<float>(current_context, d_img);
    else
      this->d_img_raw->copy(d_img, 1, 1);
  }

  if (this->d_dark != nullptr) d_img->axpy(-1.f, this->d_dark, 1, 1);

  if (this->d_flat != nullptr)
    mult_vect(d_img->getData(), this->d_flat->getData(), d_img->getNbElem(),
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

int sutra_centroider::set_npix(int npix) {
  this->npix = npix;

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

int sutra_centroider::set_centroids_ref(float *centroids_ref) {
  this->d_centroids_ref->host2device(centroids_ref);
  return EXIT_SUCCESS;
}
