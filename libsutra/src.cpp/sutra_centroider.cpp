#include <sutra_centroider.h>

template <class Tin, class Tout>
sutra_centroider<Tin, Tout>::sutra_centroider(carma_context *context,
                                              sutra_wfs *wfs, long nvalid,
                                              float offset, float scale,
                                              bool filter_TT, int device) {
  this->current_context = context;
  this->device = device;
  context->set_activeDevice(device, 1);
  this->wfs = wfs;

  this->nvalid = nvalid;
  this->offset = offset;
  this->scale = scale;
  this->nslopes = 0;
  this->npix = 0;
  this->nxsub = 0;
  this->filter_TT = filter_TT;

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
  this->d_lutPix = nullptr;
  this->d_bincube = nullptr;
  this->d_validMask = nullptr;
  this->d_centro_filtered = nullptr;
  this->d_ref_Tip = nullptr;
  this->d_ref_Tilt = nullptr;
  this->d_TT_slopes = nullptr;
}

template <class Tin, class Tout>
sutra_centroider<Tin, Tout>::~sutra_centroider() {
  if (this->d_intensities != nullptr) delete this->d_intensities;
  if (this->d_centroids_ref != nullptr) delete this->d_centroids_ref;
  if (this->d_img != nullptr) delete this->d_img;
  if (this->d_img_raw != nullptr) delete this->d_img_raw;
  if (this->d_validx != nullptr) delete this->d_validx;
  if (this->d_validy != nullptr) delete this->d_validy;
  if (this->d_dark != nullptr) delete this->d_dark;
  if (this->d_flat != nullptr) delete this->d_flat;
  if (this->d_lutPix != nullptr) delete this->d_lutPix;
  if (this->d_bincube != nullptr) delete this->d_bincube;
  if (this->d_validMask != nullptr) delete this->d_validMask;
  if (this->d_TT_slopes != nullptr) delete this->d_TT_slopes;
  if (this->d_centro_filtered != nullptr) delete this->d_centro_filtered;
  if (this->d_ref_Tip != nullptr) delete this->d_ref_Tip;
  if (this->d_ref_Tilt != nullptr) delete this->d_ref_Tilt;
}

template <class Tin, class Tout>
int sutra_centroider<Tin, Tout>::set_scale(float scale) {
  this->scale = scale;
  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int sutra_centroider<Tin, Tout>::set_offset(float offset) {
  this->offset = offset;
  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int sutra_centroider<Tin, Tout>::set_nxsub(int nxsub) {
  this->nxsub = nxsub;
  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int sutra_centroider<Tin, Tout>::init_calib(int n, int m) {
  current_context->set_activeDevice(device, 1);
  if (this->d_dark == nullptr) {
    long dims_data2[3] = {2, n, m};
    this->d_dark = new carma_obj<float>(current_context, dims_data2);
    this->d_dark->reset();
  }
  if (this->d_flat == nullptr) {
    long dims_data2[3] = {2, n, m};
    this->d_flat = new carma_obj<float>(current_context, dims_data2);
    this->d_flat->memSet(1.f);
  }
  if (this->d_lutPix == nullptr) {
    long dims_data1[3] = {1, n * m};
    this->d_lutPix = new carma_obj<int>(current_context, dims_data1);
    std::vector<int> h_lutPix(n * m);
    for (int i = 0; i < n * m; ++i)
      h_lutPix[i] = i;
    this->d_lutPix->host2device(h_lutPix.data());
  }
  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int sutra_centroider<Tin, Tout>::init_roi(int N) {
  current_context->set_activeDevice(device, 1);
  if (this->d_validx == nullptr) {
    long dims_data[2] = {1, N};
    this->d_validx = new carma_obj<int>(current_context, dims_data);
    this->d_validy = new carma_obj<int>(current_context, dims_data);
  }
  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int sutra_centroider<Tin, Tout>::set_dark(float *dark, int n) {
  current_context->set_activeDevice(device, 1);
  if (this->d_dark == nullptr) {
    long dims_data2[3] = {2, n, n};
    this->d_dark = new carma_obj<float>(current_context, dims_data2);
  }
  this->d_dark->host2device(dark);
  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int sutra_centroider<Tin, Tout>::set_flat(float *flat, int n) {
  current_context->set_activeDevice(device, 1);
  if (this->d_flat == nullptr) {
    long dims_data2[3] = {2, n, n};
    this->d_flat = new carma_obj<float>(current_context, dims_data2);
  }
  this->d_flat->host2device(flat);
  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int sutra_centroider<Tin, Tout>::set_lutPix(int *lutPix, int n) {
  current_context->set_activeDevice(device, 1);
  if (this->d_lutPix == nullptr)
  {
    long dims_data1[2] = {1, n};
    this->d_lutPix = new carma_obj<int>(current_context, dims_data1);
  }
  this->d_lutPix->host2device(lutPix);
  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int sutra_centroider<Tin, Tout>::get_validMask()
{
  this->current_context->set_activeDevice(this->device, 1);
  if (this->d_validMask == nullptr) {
    if (this->d_img == nullptr) {
      std::cout << "RTC image has not been initialized" << std::endl;
      return EXIT_FAILURE;
    }
    this->d_validMask = new carma_obj<int>(current_context, d_img->getDims());
    this->d_validMask->reset();
  }

  fill_validMask(this->d_validMask->getDims(1), this->npix, this->nvalid,
                 this->d_validMask->getData(), this->d_validx->getData(),
                 this->d_validy->getData(),
                 current_context->get_device(device));

  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int sutra_centroider<Tin, Tout>::calibrate_img() {
  current_context->set_activeDevice(device, 1);

  if (this->d_img_raw == nullptr)
  {
    std::cout << "Image not initialized\n"
              << std::endl;
    return EXIT_FAILURE;
  }

  const long *dims = this->d_img->getDims();
  init_calib(dims[0], dims[1]);

  calibration<Tin>(this->d_img_raw->getData(), this->d_img->getData(),
                   this->d_dark->getData(), this->d_flat->getData(),
                   this->d_lutPix->getData(), this->d_img->getNbElem(),
                   this->current_context->get_device(this->device));

  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int sutra_centroider<Tin, Tout>::load_img(carma_obj<Tin> *img) {
  return this->load_img(img->getData(), img->getDims(1), img->getDevice());
}

template <class Tin, class Tout>
int sutra_centroider<Tin, Tout>::load_img(Tin *img, int n) {
  return this->load_img(img, n, -1);
}

template <class Tin, class Tout>
int sutra_centroider<Tin, Tout>::load_img(Tin *img, int n, int location) {
  current_context->set_activeDevice(device, 1);
  if (this->d_img_raw == nullptr) {
    long dims_data2[3] = {2, n, n};
    this->d_img_raw = new carma_obj<Tin>(current_context, dims_data2);
    this->d_img = new carma_obj<float>(current_context, dims_data2);
  }

  if (location < 0) {  // img data on host
    this->d_img_raw->host2device(img);
  } else {  // img data on device
    this->d_img_raw->copyFrom(img, this->d_img_raw->getNbElem());
  }
  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int sutra_centroider<Tin, Tout>::set_npix(int npix) {
  this->npix = npix;

  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int sutra_centroider<Tin, Tout>::load_validpos(int *ivalid, int *jvalid,
                                               int N) {
  current_context->set_activeDevice(device, 1);
  if (this->d_validx == nullptr) {
    this->init_roi(N);
  }

  this->d_validx->host2device(ivalid);
  this->d_validy->host2device(jvalid);

  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int sutra_centroider<Tin, Tout>::set_centroids_ref(float *centroids_ref) {
  this->d_centroids_ref->host2device(centroids_ref);
  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int sutra_centroider<Tin, Tout>::init_TT_filter() {
  this->current_context->set_activeDevice(device, 1);
  long dims_data[2] = {1, 2};
  this->d_TT_slopes = new carma_obj<float>(this->current_context, dims_data);
  dims_data[1] = this->nslopes;
  this->d_centro_filtered =
      new carma_obj<float>(this->current_context, dims_data);
  this->d_ref_Tip = new carma_obj<float>(this->current_context, dims_data);
  this->d_ref_Tilt = new carma_obj<float>(this->current_context, dims_data);

  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int sutra_centroider<Tin, Tout>::apply_TT_filter(Tout *centroids) {
  return this->apply_TT_filter_impl(centroids, std::is_same<Tout, float>());
}

template <class Tin, class Tout>
template <typename Q>
typename std::enable_if<std::is_same<Q, float>::value, int>::type
sutra_centroider<Tin, Tout>::apply_TT_filter_impl(Tout *centroids,
                                                  std::true_type) {
  this->d_centro_filtered->copyFrom(centroids, this->nslopes);

  float tip = this->d_centro_filtered->dot(this->d_ref_Tip, 1, 1);
  float tilt = this->d_centro_filtered->dot(this->d_ref_Tilt, 1, 1);

  this->d_centro_filtered->axpy(-1.f * tip, this->d_ref_Tip, 1, 1);
  this->d_centro_filtered->axpy(-1.f * tilt, this->d_ref_Tilt, 1, 1);

  this->d_centro_filtered->copyInto(centroids, this->nslopes);

  float TT_data[2] = {tip, tilt};
  this->d_TT_slopes->host2device(TT_data);

  return EXIT_SUCCESS;
}

template <class Tin, class Tout>
int sutra_centroider<Tin, Tout>::apply_TT_filter_impl(Tout *centroids,
                                                      std::false_type) {
  DEBUG_TRACE("Tip/tilt filtering is only implemented in single precision");
  return EXIT_SUCCESS;
}

template class sutra_centroider<float, float>;
template class sutra_centroider<uint16_t, float>;
#ifdef CAN_DO_HALF
template class sutra_centroider<float, half>;
template class sutra_centroider<uint16_t, half>;
#endif
