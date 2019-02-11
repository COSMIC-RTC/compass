#include <cured.h>
#include <sutra_controller_cured.h>
#include <string>

template <typename T, typename Tout>
sutra_controller_cured<T, Tout>::sutra_controller_cured(
    carma_context *context, long nvalid, long nslopes, long nactu, float delay,
    sutra_dms *dms, int *idx_dms, int ndm)
    : sutra_controller<T, Tout>(context, nvalid, nslopes, nactu, delay, dms,
                                idx_dms, ndm),
      gain(0),
      ndivs(0),
      tt_flag(false),
      h_syscure(nullptr),
      h_parcure(nullptr) {
  long dims_data1[2] = {1, 0};
  long dims_data2[3] = {2, 0, 0};

  dims_data1[1] = nslopes;
  this->d_centroids = new carma_obj<T>(context, dims_data1);

  dims_data1[1] = nslopes;
  this->h_centroids = new carma_host_obj<T>(dims_data1, MA_PAGELOCK);

  dims_data1[1] = nactu;
  this->h_err = new carma_host_obj<T>(dims_data1, MA_PAGELOCK);
  this->d_err = new carma_obj<T>(context, dims_data1);

  dims_data2[1] = nslopes;
  dims_data2[2] = nactu;

  this->d_imat = new carma_obj<T>(context, dims_data2);
  if ((int)delay > 0) {
    dims_data2[1] = nslopes;
    dims_data2[2] = (int)delay + 1;
    this->d_cenbuff = new carma_obj<T>(context, dims_data2);
  } else
    this->d_cenbuff = 0L;
}

template <typename T, typename Tout>
sutra_controller_cured<T, Tout>::~sutra_controller_cured() {
  this->current_context->set_activeDevice(this->device, 1);

  if (this->h_centroids != nullptr) delete this->h_centroids;
  if (this->h_err != nullptr) delete this->h_err;
  if (this->d_imat != nullptr) delete this->d_imat;
  if (this->d_err != nullptr) delete this->d_err;
  if (this->d_cenbuff) delete this->d_cenbuff;
  curefree((sysCure *)this->h_syscure, (parCure *)this->h_parcure);
}

template <typename T, typename Tout>
int sutra_controller_cured<T, Tout>::set_gain(T gain) {
  this->gain = gain;
  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int sutra_controller_cured<T, Tout>::comp_com() {
  this->current_context->set_activeDevice(this->device, 1);

  // this->frame_delay();
  h_centroids->cpy_obj(this->d_centroids, cudaMemcpyDeviceToHost);

  if (this->tt_flag) {
    cured((sysCure *)this->h_syscure, (parCure *)this->h_parcure,
          this->h_centroids->getData(), this->h_err->getData(),
          this->h_err->getDataAt(this->h_err->getNbElem() - 2),
          this->h_err->getDataAt(this->h_err->getNbElem() - 1));

    //*this->h_err->getDataAt(this->h_err->getNbElem()-2) *= -1.0f;
    //*this->h_err->getDataAt(this->h_err->getNbElem()-1) *= -1.0f;
  } else {
    cured((sysCure *)this->h_syscure, (parCure *)this->h_parcure,
          this->h_centroids->getData(), this->h_err->getData());
  }

  h_err->cpy_obj(this->d_err, cudaMemcpyHostToDevice);

  mult_int(this->d_com->getData(), this->d_err->getData(), -1.0f * this->gain,
           this->nactu(), this->current_context->get_device(this->device));

  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int sutra_controller_cured<T, Tout>::init_cured(int nxsubs, int *isvalid,
                                                int ndivs, int tt) {
  if (tt > 0)
    this->tt_flag = true;
  else
    this->tt_flag = false;
  this->ndivs = (ndivs > 0) ? ndivs : 1;
  this->h_syscure = (void *)cureSystem(
      nxsubs, this->nslope() / 2.,
      this->tt_flag ? this->nactu() - 2 : this->nactu(), isvalid, this->ndivs);
  this->h_parcure = (void *)cureInit((sysCure *)this->h_syscure);
  return EXIT_SUCCESS;
}

template <typename T, typename Tout>
int sutra_controller_cured<T, Tout>::frame_delay() {
  // here we place the content of d_centroids into cenbuf and get
  // the actual centroid frame for error computation depending on delay value

  if ((int)this->delay > 0) {
    for (int cc = 0; cc < this->delay; cc++)
      shift_buf(&((this->d_cenbuff->getData())[cc * this->nslope()]), 1,
                this->nslope(),
                this->current_context->get_device(this->device));

    carmaSafeCall(cudaMemcpy(
        this->d_cenbuff->getDataAt((int)this->delay * this->nslope()),
        this->d_centroids->getData(), sizeof(T) * this->nslope(),
        cudaMemcpyDeviceToDevice));

    carmaSafeCall(
        cudaMemcpy(this->d_centroids->getData(), this->d_cenbuff->getData(),
                   sizeof(T) * this->nslope(), cudaMemcpyDeviceToDevice));
  }

  return EXIT_SUCCESS;
}

template class sutra_controller_cured<float, float>;
template class sutra_controller_cured<float, uint16_t>;
