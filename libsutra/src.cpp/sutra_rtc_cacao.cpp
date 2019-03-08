#ifndef USE_OCTOPUS

#include <sutra_rtc_cacao.h>

template <typename Tin, typename Tcomp, typename Tout>
sutra_rtc_cacao<Tin, Tcomp, Tout>::sutra_rtc_cacao(std::string iCalFrame_name,
                                                   std::string iLoopFrame_name)
    : sutra_rtc<Tin, Tcomp, Tout>(),
      iCalFrame_name_(iCalFrame_name),
      iLoopFrame_name_(iLoopFrame_name),
      framecounter_(0),
      nslp_(0),
      ncmd_(0),
      nvalid_(0),
      is_initialised_(false) {
  iCalFrame_.reset();
  iLoopFrame_.reset();
}

template <typename Tin, typename Tcomp, typename Tout>
sutra_rtc_cacao<Tin, Tcomp, Tout>::~sutra_rtc_cacao() {
  if (!is_initialised_) {
    return;
  }
  iCalFrame_.reset();
  iLoopFrame_.reset();
}

template <typename Tin, typename Tcomp, typename Tout>
void sutra_rtc_cacao<Tin, Tcomp, Tout>::allocateBuffers() {
  if (is_initialised_) {
    return;
  }

  uint32_t size = this->d_centro[0]->d_img->getDims(1);
  iCalFrame_ = std::make_shared<ipc::Cacao<float>>(ipc::Cacao<float>(
      iCalFrame_name_, std::vector<uint32_t>{10, size, size}));

  nslp_ = 0;
  ncmd_ = 0;
  for (unsigned int i = 0; i < this->d_control.size(); i++) {
    nslp_ += this->d_control[i]->nslope();
    ncmd_ += this->d_control[i]->nactu();
  }
  nvalid_ = nslp_ / 2;

  uint32_t size_tot = nvalid_ + nslp_ + ncmd_;
  iLoopFrame_ = std::make_shared<ipc::Cacao<Tcomp>>(ipc::Cacao<Tcomp>(
      iLoopFrame_name_, std::vector<uint32_t>{10, 1, size_tot}));

  is_initialised_ = true;
}

template <typename Tin, typename Tcomp, typename Tout>
void sutra_rtc_cacao<Tin, Tcomp, Tout>::publish() {
  if (!is_initialised_) {
    allocateBuffers();
  }

  float* iFrame = iCalFrame_->outputPtr();

  this->d_centro[0]->d_img->device2host(iFrame);
  iCalFrame_->notify();

  Tcomp* zFrame = iLoopFrame_->outputPtr();

  for (unsigned int i = 0; i < this->d_centro.size(); i++) {
    this->d_centro[i]->d_intensities->device2host(zFrame);
    zFrame += this->d_centro[i]->nvalid;
  }

  for (unsigned int i = 0; i < this->d_control.size(); i++) {
    this->d_control[i]->d_centroids->device2host(zFrame);
    zFrame += this->d_control[i]->nslope();
  }

  for (unsigned int i = 0; i < this->d_control.size(); i++) {
    this->d_control[i]->d_com->device2host(zFrame);
    zFrame += this->d_control[i]->nactu();
  }

  iLoopFrame_->notify();

  framecounter_++;
}
template class sutra_rtc_cacao<float, float, float>;
template class sutra_rtc_cacao<uint16_t, float, float>;
template class sutra_rtc_cacao<float, float, uint16_t>;
template class sutra_rtc_cacao<uint16_t, float, uint16_t>;
#ifdef CAN_DO_HALF
// template class sutra_rtc_cacao<half>;
#endif

#endif /* USE_BRAHMA */
