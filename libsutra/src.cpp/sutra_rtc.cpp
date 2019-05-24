#include <sutra_rtc.h>
#include <sutra_wfs_pyr_pyrhr.h>
#include <sutra_wfs_sh.h>

template <typename Tin, typename T, typename Tout>
sutra_rtc<Tin, T, Tout>::sutra_rtc() {}

template <typename Tin, typename T, typename Tout>
sutra_rtc<Tin, T, Tout>::~sutra_rtc() {
  //  for (size_t idx = 0; idx < (this->d_centro).size(); idx++) {
  while ((this->d_centro).size() > 0) {
    delete this->d_centro.back();
    (this->d_centro).pop_back();
  }

  //  for (size_t idx = 0; idx < (this->d_control).size(); idx++) {
  while ((this->d_control).size() > 0) {
    delete this->d_control.back();
    (this->d_control).pop_back();
  }

  // delete this->context;
}

template <typename Tin, typename T, typename Tout>
int sutra_rtc<Tin, T, Tout>::remove_centroider(int ncentro) {
  delete this->d_centro[ncentro];
  this->d_centro.erase(this->d_centro.begin() + ncentro);
  return EXIT_SUCCESS;
}

template <typename Tin, typename T, typename Tout>
int sutra_rtc<Tin, T, Tout>::remove_controller(int ncontrol) {
  delete this->d_control[ncontrol];
  this->d_control.erase(this->d_control.begin() + ncontrol);
  return EXIT_SUCCESS;
}

template <typename Tin, typename T, typename Tout>
int sutra_rtc<Tin, T, Tout>::add_centroider(carma_context *context, long nvalid,
                                            float offset, float scale,
                                            bool filter_TT, long device,
                                            std::string typec) {
  return add_centroider(context, nvalid, offset, scale, filter_TT, device,
                        typec, nullptr);
}

template <typename Tin, typename T, typename Tout>
int sutra_rtc<Tin, T, Tout>::add_centroider(carma_context *context, long nvalid,
                                            float offset, float scale,
                                            bool filter_TT, long device,
                                            std::string typec, sutra_wfs *wfs) {
  return add_centroider_impl(context, nvalid, offset, scale, filter_TT, device,
                             typec, wfs, std::is_same<T, half>());
}

template <typename Tin, typename T, typename Tout>
template <typename Q>
typename std::enable_if<!std::is_same<Q, half>::value, int>::type
sutra_rtc<Tin, T, Tout>::add_centroider_impl(carma_context *context,
                                             long nvalid, float offset,
                                             float scale, bool filter_TT,
                                             long device, std::string typec,
                                             sutra_wfs *wfs, std::false_type) {
  if (typec.compare("bpcog") == 0)
    this->d_centro.push_back(new sutra_centroider_bpcog<Tin, T>(
        context, wfs, nvalid, offset, scale, filter_TT, device, 10));
  else if (typec.compare("cog") == 0)
    this->d_centro.push_back(new sutra_centroider_cog<Tin, T>(
        context, wfs, nvalid, offset, scale, filter_TT, device));
  else if (typec.compare("corr") == 0)
    this->d_centro.push_back(new sutra_centroider_corr<Tin, T>(
        context, wfs, nvalid, offset, scale, filter_TT, device));
  else if (typec.compare("pyr") == 0)
    this->d_centro.push_back(new sutra_centroider_pyr<Tin, T>(
        context, wfs, nvalid, offset, scale, filter_TT, device));
  else if (typec.compare("tcog") == 0)
    this->d_centro.push_back(new sutra_centroider_tcog<Tin, T>(
        context, wfs, nvalid, offset, scale, filter_TT, device));
  else if (typec.compare("wcog") == 0)
    this->d_centro.push_back(new sutra_centroider_wcog<Tin, T>(
        context, wfs, nvalid, offset, scale, filter_TT, device));
  else if (typec.compare("maskedpix") == 0) {
    if (wfs == nullptr) {
      this->d_centro.push_back(new sutra_centroider_maskedPix<Tin, T>(
          context, wfs, nvalid, 4, offset, scale, filter_TT, device));
    } else if (wfs->type == "pyrhr") {
      sutra_wfs_pyr_pyrhr *pwfs = dynamic_cast<sutra_wfs_pyr_pyrhr *>(wfs);
      this->d_centro.push_back(new sutra_centroider_maskedPix<Tin, T>(
          context, pwfs, nvalid, pwfs->npupils, offset, scale, filter_TT,
          device));
    } else
      DEBUG_TRACE("WFS must be pyrhr");
  } else {
    std::cerr << "centroider unknown" << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

#ifdef CAN_DO_HALF
template <typename Tin, typename T, typename Tout>
int sutra_rtc<Tin, T, Tout>::add_centroider_impl(carma_context *context,
                                                 long nvalid, float offset,
                                                 float scale, bool filter_TT,
                                                 long device, std::string typec,
                                                 sutra_wfs *wfs,
                                                 std::true_type) {
  if (typec.compare("cog") == 0)
    this->d_centro.push_back(new sutra_centroider_cog<Tin, T>(
        context, wfs, nvalid, offset, scale, filter_TT, device));
  else if (typec.compare("bpcog") == 0)
    this->d_centro.push_back(new sutra_centroider_bpcog<Tin, T>(
        context, wfs, nvalid, offset, scale, filter_TT, device, 10));
  else if (typec.compare("tcog") == 0)
    this->d_centro.push_back(new sutra_centroider_tcog<Tin, T>(
        context, wfs, nvalid, offset, scale, filter_TT, device));
  else if (typec.compare("maskedpix") == 0) {
    if (wfs == nullptr) {
      this->d_centro.push_back(new sutra_centroider_maskedPix<Tin, T>(
          context, wfs, nvalid, 4, offset, scale, filter_TT, device));
    } else if (wfs->type == "pyrhr") {
      sutra_wfs_pyr_pyrhr *pwfs = dynamic_cast<sutra_wfs_pyr_pyrhr *>(wfs);
      this->d_centro.push_back(new sutra_centroider_maskedPix<Tin, T>(
          context, pwfs, nvalid, pwfs->npupils, offset, scale, filter_TT,
          device));
    } else
      DEBUG_TRACE("WFS must be pyrhr");
  } else
    DEBUG_TRACE("Not implemented for half precision yet");

  return EXIT_SUCCESS;
}
#endif

template <typename Tin, typename T, typename Tout>
int sutra_rtc<Tin, T, Tout>::add_controller(carma_context *context, int nvalid,
                                            int nslope, int nactu, float delay,
                                            long device, std::string typec,
                                            sutra_dms *dms, int *idx_dms,
                                            int ndm, int Nphi,
                                            bool wfs_direction) {
  return add_controller_impl(context, this->d_control, nvalid, nslope, nactu,
                             delay, device, typec, dms, idx_dms, ndm, Nphi,
                             wfs_direction, std::is_same<T, half>());
}

template <typename Tin, typename T, typename Tout>
template <typename Q>
typename std::enable_if<!std::is_same<Q, half>::value, int>::type
sutra_rtc<Tin, T, Tout>::add_controller_impl(
    carma_context *context, vector<sutra_controller<T, Tout> *> &d_control,
    int nvalid, int nslope, int nactu, float delay, long device,
    std::string typec, sutra_dms *dms, int *idx_dms, int ndm, int Nphi,
    bool wfs_direction, std::false_type) {
  if (typec.compare("ls") == 0) {
    d_control.push_back(new sutra_controller_ls<T, Tout>(
        context, nvalid, nslope, nactu, delay, dms, idx_dms, ndm));
  } else if (typec.compare("geo") == 0) {
    d_control.push_back(new sutra_controller_geo<T, Tout>(
        context, nactu, Nphi, delay, dms, idx_dms, ndm, wfs_direction));

  } else if (typec.compare("cured") == 0) {
    d_control.push_back(new sutra_controller_cured<T, Tout>(
        context, nvalid, nslope, nactu, delay, dms, idx_dms, ndm));
  } else if (typec.compare("mv") == 0) {
    d_control.push_back(new sutra_controller_mv<T, Tout>(
        context, nvalid, nslope, nactu, delay, dms, idx_dms, ndm));
  } else if (typec.compare("generic") == 0) {
    d_control.push_back(new sutra_controller_generic<T, Tout>(
        context, nvalid, nslope, nactu, delay, dms, idx_dms, ndm));
    // } else if ((typec.compare("kalman_GPU") == 0) ||
    //            (typec.compare("kalman_CPU") == 0)) {
    //   d_control.push_back(
    //       new sutra_controller_kalman(context, nslope, nactu, dms, idx_dms,
    //       ndm));
  } else {
    DEBUG_TRACE("Controller '%s' unknown\n", typec.c_str());
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

template <typename Tin, typename T, typename Tout>
int sutra_rtc<Tin, T, Tout>::add_controller_impl(
    carma_context *context, vector<sutra_controller<T, Tout> *> &d_control,
    int nvalid, int nslope, int nactu, float delay, long device,
    std::string typec, sutra_dms *dms, int *idx_dms, int ndm, int Nphi,
    bool wfs_direction, std::true_type) {
  if (typec.compare("generic") == 0) {
    d_control.push_back(new sutra_controller_generic<T, Tout>(
        context, nvalid, nslope, nactu, delay, dms, idx_dms, ndm));

  } else {
    DEBUG_TRACE("Not implemented in half precision yet");
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

template <typename Tin, typename T, typename Tout>
int sutra_rtc<Tin, T, Tout>::do_imat(int ncntrl, sutra_dms *ydm) {
  return do_imat_impl(ncntrl, ydm, std::is_same<T, float>());
}

template <typename Tin, typename T, typename Tout>
int sutra_rtc<Tin, T, Tout>::do_imat_impl(int ncntrl, sutra_dms *ydm,
                                          std::false_type) {
  DEBUG_TRACE("Not implemented for this computation type");
  return EXIT_FAILURE;
}

template <typename Tin, typename T, typename Tout>
template <typename Q>
typename std::enable_if<std::is_same<Q, float>::value, int>::type
sutra_rtc<Tin, T, Tout>::do_imat_impl(int ncntrl, sutra_dms *ydm,
                                      std::true_type) {
  carma_obj<T> *d_imat = NULL;
  if (this->d_control[ncntrl]->get_type().compare("ls") == 0) {
    sutra_controller_ls<T, Tout> *control =
        dynamic_cast<sutra_controller_ls<T, Tout> *>(this->d_control[ncntrl]);
    d_imat = control->d_imat;
  } else if (this->d_control[ncntrl]->get_type().compare("mv") == 0) {
    sutra_controller_mv<T, Tout> *control =
        dynamic_cast<sutra_controller_mv<T, Tout> *>(this->d_control[ncntrl]);
    d_imat = control->d_imat;
  } else {
    DEBUG_TRACE("Controller needs to be LS or MV\n");
    return EXIT_SUCCESS;
  }

  vector<sutra_dm *>::iterator p;

  p = this->d_control[ncntrl]->d_dmseen.begin();

  while (p != this->d_control[ncntrl]->d_dmseen.end()) {
    sutra_dm *dm = *p;
    if (dm->type == "tt") {
      for (size_t idx_cntr = 0; idx_cntr < (this->d_centro).size();
           idx_cntr++) {
        if (this->d_centro[idx_cntr]->filter_TT) {
          std::cout << "Measuring TT reference for centro : " << idx_cntr
                    << std::endl;
          this->d_centro[idx_cntr]->filter_TT = false;

          // Tip Push
          dm->comp_oneactu(0, dm->push4imat);

          this->comp_images_imat(ydm);
          this->d_centro[idx_cntr]->get_cog(
              this->d_centro[idx_cntr]->d_intensities->getData(),
              this->d_centro[idx_cntr]->d_centro_filtered->getData(), true);
          this->d_centro[idx_cntr]->d_centro_filtered->scale(
              0.5f / dm->push4imat, 1);
          this->d_centro[idx_cntr]->d_centro_filtered->copyInto(
              this->d_centro[idx_cntr]->d_ref_Tip->getData(),
              this->d_centro[idx_cntr]->nslopes);

          dm->reset_shape();

          // Tip Pull
          dm->comp_oneactu(0, -1.0f * dm->push4imat);
          this->comp_images_imat(ydm);
          this->d_centro[idx_cntr]->get_cog(
              this->d_centro[idx_cntr]->d_intensities->getData(),
              this->d_centro[idx_cntr]->d_centro_filtered->getData(), true);
          float alphai = -.5f / dm->push4imat;

          this->d_centro[idx_cntr]->d_ref_Tip->axpy(
              T(alphai), this->d_centro[idx_cntr]->d_centro_filtered, 1, 1);

          dm->reset_shape();

          // Tilt Push
          dm->comp_oneactu(1, dm->push4imat);

          this->comp_images_imat(ydm);
          this->d_centro[idx_cntr]->get_cog(
              this->d_centro[idx_cntr]->d_intensities->getData(),
              this->d_centro[idx_cntr]->d_centro_filtered->getData(), true);
          this->d_centro[idx_cntr]->d_centro_filtered->scale(
              0.5f / dm->push4imat, 1);
          this->d_centro[idx_cntr]->d_centro_filtered->copyInto(
              this->d_centro[idx_cntr]->d_ref_Tilt->getData(),
              this->d_centro[idx_cntr]->nslopes);

          dm->reset_shape();

          // Tilt Pull
          dm->comp_oneactu(1, -1.0f * dm->push4imat);
          this->comp_images_imat(ydm);
          this->d_centro[idx_cntr]->get_cog(
              this->d_centro[idx_cntr]->d_intensities->getData(),
              this->d_centro[idx_cntr]->d_centro_filtered->getData(), true);

          this->d_centro[idx_cntr]->d_ref_Tilt->axpy(
              T(alphai), this->d_centro[idx_cntr]->d_centro_filtered, 1, 1);

          dm->reset_shape();

          this->d_centro[idx_cntr]->filter_TT = true;

          float nrmTip = this->d_centro[idx_cntr]->d_ref_Tip->nrm2(1);
          float nrmTilt = this->d_centro[idx_cntr]->d_ref_Tilt->nrm2(1);

          this->d_centro[idx_cntr]->d_ref_Tip->scale(T(1.0f / nrmTip), 1);
          this->d_centro[idx_cntr]->d_ref_Tilt->scale(T(1.0f / nrmTilt), 1);
        }
      }
    }
    ++p;
  }

  p = this->d_control[ncntrl]->d_dmseen.begin();
  int inds1 = 0;
  int cc2 = 0;

  std::cout << "Doing imat..." << std::endl;
  while (p != this->d_control[ncntrl]->d_dmseen.end()) {
    sutra_dm *dm = *p;
    std::string desc = "DM" + carma_utils::to_string(cc2);
    auto progress = carma_utils::ProgressBar(dm->nactus, desc);
    for (int j = 0; j < dm->nactus; ++j) {
      // Push
      dm->comp_oneactu(j, dm->push4imat);

      this->comp_images_imat(ydm);
      do_centroids(ncntrl, true);

      int device = this->d_control[ncntrl]->d_centroids->getDevice();
      this->d_control[ncntrl]->d_centroids->scale(0.5f / dm->push4imat, 1);
      this->d_control[ncntrl]->d_centroids->copyInto(
          d_imat->getDataAt(inds1), this->d_control[ncntrl]->nslope());
      dm->reset_shape();
      // Pull
      dm->comp_oneactu(j, -1.0f * dm->push4imat);
      this->comp_images_imat(ydm);
      device = this->d_control[ncntrl]->d_centroids->getDevice();
      do_centroids(ncntrl, true);

      float alphai = -.5f / dm->push4imat;
      cublasSaxpy(this->d_control[ncntrl]->current_context->get_cublasHandle(),
                  this->d_control[ncntrl]->d_centroids->getNbElem(), &alphai,
                  this->d_control[ncntrl]->d_centroids->getData(), 1,
                  d_imat->getDataAt(inds1), 1);

      dm->reset_shape();
      inds1 += this->d_control[ncntrl]->nslope();
      progress.update();
      // printf("\rDoing imat...%d%%",(cc*100/nactu));
    }
    progress.finish();
    ++p;
    ++cc2;
  }
  printf("\n");
  return EXIT_SUCCESS;
}

template <typename Tin, typename T, typename Tout>
int sutra_rtc<Tin, T, Tout>::do_imat_basis(int ncntrl, sutra_dms *ydm,
                                           int nModes, T *m2v, T *pushAmpl) {
  return do_imat_basis_impl(ncntrl, ydm, nModes, m2v, pushAmpl,
                            std::is_same<T, float>());
}

template <typename Tin, typename T, typename Tout>
int sutra_rtc<Tin, T, Tout>::do_imat_basis_impl(int ncntrl, sutra_dms *ydm,
                                                int nModes, T *m2v, T *pushAmpl,
                                                std::false_type) {
  DEBUG_TRACE("Not implemented for this computation type");
  return EXIT_FAILURE;
}

template <typename Tin, typename T, typename Tout>
template <typename Q>
typename std::enable_if<std::is_same<Q, float>::value, int>::type
sutra_rtc<Tin, T, Tout>::do_imat_basis_impl(int ncntrl, sutra_dms *ydm,
                                            int nModes, T *m2v, T *pushAmpl,
                                            std::true_type) {
  carma_obj<T> d_m2v(this->d_control[ncntrl]->current_context,
                     std::vector<long>{2, ydm->nact_total(), nModes}.data(),
                     m2v);
  carma_obj<T> d_comm(this->d_control[ncntrl]->current_context,
                      std::vector<long>{1, ydm->nact_total()}.data());

  carma_obj<T> *d_imat = NULL;

  if (this->d_control[ncntrl]->get_type().compare("ls") == 0) {
    sutra_controller_ls<T, Tout> *control =
        dynamic_cast<sutra_controller_ls<T, Tout> *>(this->d_control[ncntrl]);
    long dims[3] = {2, control->d_imat->getDims()[1], nModes};
    d_imat = new carma_obj<T>(this->d_control[ncntrl]->current_context, dims);
    delete control->d_imat;
    control->d_imat = d_imat;
  } else if (this->d_control[ncntrl]->get_type().compare("mv") == 0) {
    sutra_controller_mv<T, Tout> *control =
        dynamic_cast<sutra_controller_mv<T, Tout> *>(this->d_control[ncntrl]);
    long dims[3] = {2, control->d_imat->getDims()[1], nModes};
    d_imat = new carma_obj<T>(this->d_control[ncntrl]->current_context, dims);
    delete control->d_imat;
    control->d_imat = d_imat;
  } else {
    DEBUG_TRACE("Controller needs to be LS or MV\n");
    return EXIT_SUCCESS;
  }
  int inds1 = 0;

  std::cout << "Doing imat modal..." << std::endl;
  std::string desc = "Modal iMat";
  auto progress = carma_utils::ProgressBar(nModes, desc);
  for (int j = 0; j < nModes; ++j) {
    // For each mode
    d_comm.copyFrom(d_m2v.getDataAt(j * ydm->nact_total()), ydm->nact_total());
    d_comm.scale(pushAmpl[j], 1);
    // Push
    int actuCount = 0;
    vector<sutra_dm *>::iterator p = this->d_control[ncntrl]->d_dmseen.begin();
    while (p != this->d_control[ncntrl]->d_dmseen.end()) {
      // Set each dm
      sutra_dm *dm = *p;
      dm->comp_shape(d_comm.getDataAt(actuCount));
      actuCount += dm->nactus;
      ++p;
    }
    this->comp_images_imat(ydm);  // Raytrace & compute all WFS
    do_centroids(ncntrl, true);
    int device = this->d_control[ncntrl]->d_centroids->getDevice();
    this->d_control[ncntrl]->d_centroids->scale(0.5f / pushAmpl[j], 1);
    this->d_control[ncntrl]->d_centroids->copyInto(
        d_imat->getDataAt(inds1), this->d_control[ncntrl]->nslope());

    // Pull
    actuCount = 0;
    d_comm.scale(-1.0f, 1);
    p = this->d_control[ncntrl]->d_dmseen.begin();
    while (p != this->d_control[ncntrl]->d_dmseen.end()) {
      // Set each dm
      sutra_dm *dm = *p;
      dm->comp_shape(d_comm.getDataAt(actuCount));
      actuCount += dm->nactus;
      ++p;
    }
    this->comp_images_imat(ydm);  // Raytrace & compute all WFS
    do_centroids(ncntrl, true);

    device = this->d_control[ncntrl]->d_centroids->getDevice();
    float alphai = -0.5f / pushAmpl[j];
    cublasSaxpy(this->d_control[ncntrl]->current_context->get_cublasHandle(),
                this->d_control[ncntrl]->d_centroids->getNbElem(), &alphai,
                this->d_control[ncntrl]->d_centroids->getData(), 1,
                d_imat->getDataAt(inds1), 1);

    inds1 += this->d_control[ncntrl]->nslope();
    progress.update();
  }
  progress.finish();
  printf("\n");
  return EXIT_SUCCESS;
}

template <typename Tin, typename T, typename Tout>
int sutra_rtc<Tin, T, Tout>::comp_images_imat(sutra_dms *ydm) {
  for (size_t idx_cntr = 0; idx_cntr < (this->d_centro).size(); idx_cntr++) {
    sutra_wfs *wfs = this->d_centro[idx_cntr]->wfs;
    float tmp_noise = wfs->noise;
    float tmp_nphot = wfs->nphot;
    wfs->nphot = wfs->nphot4imat;
    wfs->noise = -1;
    wfs->kernconv = true;
    wfs->sensor_trace(ydm, 1);
    wfs->comp_image();
    wfs->noise = tmp_noise;
    wfs->nphot = tmp_nphot;
    wfs->kernconv = false;
  }
  return EXIT_SUCCESS;
}

template <typename Tin, typename T, typename Tout>
int sutra_rtc<Tin, T, Tout>::do_imat_geom(int ncntrl, sutra_dms *ydm,
                                          int type) {
  return do_imat_geom_impl(ncntrl, ydm, type, std::is_same<T, float>());
}

template <typename Tin, typename T, typename Tout>
int sutra_rtc<Tin, T, Tout>::do_imat_geom_impl(int ncntrl, sutra_dms *ydm,
                                               int type, std::false_type) {
  DEBUG_TRACE("Not implemented for this computation type");
  return EXIT_FAILURE;
}

template <typename Tin, typename T, typename Tout>
template <typename Q>
typename std::enable_if<std::is_same<Q, float>::value, int>::type
sutra_rtc<Tin, T, Tout>::do_imat_geom_impl(int ncntrl, sutra_dms *ydm, int type,
                                           std::true_type) {
  vector<sutra_dm *>::iterator p;
  p = this->d_control[ncntrl]->d_dmseen.begin();
  int inds1, inds2;
  inds1 = 0;
  while (p != this->d_control[ncntrl]->d_dmseen.end()) {
    sutra_dm *dm = *p;
    for (int j = 0; j < dm->nactus; j++) {
      dm->comp_oneactu(j, dm->push4imat);  //
      inds2 = 0;
      for (size_t idx_cntr = 0; idx_cntr < (this->d_control).size();
           idx_cntr++) {
        carma_obj<T> *d_imat;
        if (this->d_control[ncntrl]->get_type().compare("ls") == 0) {
          sutra_controller_ls<T, Tout> *control =
              dynamic_cast<sutra_controller_ls<T, Tout> *>(
                  this->d_control[ncntrl]);
          d_imat = control->d_imat;
        } else if (this->d_control[ncntrl]->get_type().compare("mv") == 0) {
          sutra_controller_mv<T, Tout> *control =
              dynamic_cast<sutra_controller_mv<T, Tout> *>(
                  this->d_control[ncntrl]);
          d_imat = control->d_imat;
        } else {
          DEBUG_TRACE("controller should be ls or mv");
          return EXIT_FAILURE;
        }

        sutra_wfs *wfs = this->d_centro[idx_cntr]->wfs;

        wfs->sensor_trace(ydm, 1);
        // sensors->d_wfs[nwfs]->comp_image();
        if (wfs->type != "sh") {
          sutra_wfs_sh *_wfs = dynamic_cast<sutra_wfs_sh *>(wfs);
          _wfs->slopes_geom(d_imat[inds1 + inds2], type);
        } else if (wfs->type != "pyrhr") {
          sutra_wfs_pyr_pyrhr *_wfs = dynamic_cast<sutra_wfs_pyr_pyrhr *>(wfs);
          _wfs->slopes_geom(d_imat[inds1 + inds2], type);
        } else {
          DEBUG_TRACE("wfs should be a SH, a geo or a pyrhr");
          return EXIT_FAILURE;
        }
        inds2 += 2 * wfs->nvalid;
      }
      dm->reset_shape();
      inds1 += this->d_control[ncntrl]->nslope();
    }
    ++p;
  }
  return EXIT_SUCCESS;
}

template <typename Tin, typename T, typename Tout>
int sutra_rtc<Tin, T, Tout>::do_calibrate_img() {
  for (size_t idx_cntr = 0; idx_cntr < (this->d_centro).size(); idx_cntr++) {
    this->do_calibrate_img(idx_cntr);
  }

  return EXIT_SUCCESS;
}

template <typename Tin, typename T, typename Tout>
int sutra_rtc<Tin, T, Tout>::do_calibrate_img(int ncntrl) {
  this->d_centro[ncntrl]->calibrate_img();

  return EXIT_SUCCESS;
}

template <typename Tin, typename T, typename Tout>
int sutra_rtc<Tin, T, Tout>::do_centroids() {
  for (size_t idx_cntr = 0; idx_cntr < (this->d_centro).size(); idx_cntr++) {
    this->d_centro[idx_cntr]->get_cog();
  }

  return EXIT_SUCCESS;
}

template <typename Tin, typename T, typename Tout>
int sutra_rtc<Tin, T, Tout>::do_centroids(int ncntrl) {
  return do_centroids(ncntrl, true);
}

template <typename Tin, typename T, typename Tout>
int sutra_rtc<Tin, T, Tout>::do_centroids(int ncntrl, bool noise) {
  int indslope = 0;

  for (size_t idx_cntr = 0; idx_cntr < (this->d_centro).size(); idx_cntr++) {
    if (this->d_centro[idx_cntr]->wfs != nullptr) {
      this->d_centro[idx_cntr]->get_cog(
          this->d_centro[idx_cntr]->d_intensities->getData(),
          this->d_control[ncntrl]->d_centroids->getDataAt(indslope), noise);

      indslope += this->d_centro[idx_cntr]->nslopes;
    } else {
      this->d_centro[idx_cntr]->get_cog(
          this->d_centro[idx_cntr]->d_img->getData(),
          this->d_centro[idx_cntr]->d_intensities->getData(),
          this->d_control[ncntrl]->d_centroids->getDataAt(indslope),
          this->d_centro[idx_cntr]->nvalid, this->d_centro[idx_cntr]->npix,
          this->d_centro[idx_cntr]->d_img->getDims(1));

      indslope += this->d_centro[idx_cntr]->nslopes;
    }
  }
  // remove_ref(ncntrl);

  return EXIT_SUCCESS;
}

template <typename Tin, typename T, typename Tout>
int sutra_rtc<Tin, T, Tout>::do_centroids_geom(int ncntrl) {
  return do_centroids_geom_impl(ncntrl, std::is_same<T, float>());
}
template <typename Tin, typename T, typename Tout>
int sutra_rtc<Tin, T, Tout>::do_centroids_geom_impl(int ncntrl,
                                                    std::false_type) {
  DEBUG_TRACE("Not implemented for this compuation type");
  return EXIT_FAILURE;
}

template <typename Tin, typename T, typename Tout>
template <typename Q>
typename std::enable_if<std::is_same<Q, float>::value, int>::type
sutra_rtc<Tin, T, Tout>::do_centroids_geom_impl(int ncntrl, std::true_type) {
  int inds2 = 0;

  for (size_t idx_cntr = 0; idx_cntr < (this->d_centro).size(); idx_cntr++) {
    sutra_wfs *wfs = this->d_centro[idx_cntr]->wfs;
    if (wfs->type == "sh") {
      sutra_wfs_sh *_wfs = dynamic_cast<sutra_wfs_sh *>(wfs);
      _wfs->slopes_geom(this->d_control[ncntrl]->d_centroids->getDataAt(inds2));
    } else if (wfs->type == "pyrhr") {
      sutra_wfs_pyr_pyrhr *_wfs = dynamic_cast<sutra_wfs_pyr_pyrhr *>(wfs);
      _wfs->slopes_geom(this->d_control[ncntrl]->d_centroids->getDataAt(inds2));
    } else {
      DEBUG_TRACE("wfs could be a SH, geo or pyrhr");
      return EXIT_FAILURE;
    }
    /*
     this->d_centro[idx_cntr]->get_cog(sensors->d_wfs[nwfs],
     (*this->d_control[ncntrl]->d_centroids)[inds2]);
     */
    inds2 += 2 * wfs->nvalid;
  }

  return EXIT_SUCCESS;
}

template <typename Tin, typename T, typename Tout>
int sutra_rtc<Tin, T, Tout>::do_centroids_ref(int ncntrl) {
  this->do_centroids(ncntrl);
  typename vector<sutra_centroider<Tin, T> *>::iterator sc;
  sc = this->d_centro.begin();
  int inds;
  inds = 0;
  while (sc != this->d_centro.end()) {
    (*sc)->d_centroids_ref->reset();
    (*sc)->d_centroids_ref->axpy(1.0f, this->d_control[ncntrl]->d_centroids, 1,
                                 1, inds);
    inds += (*sc)->nslopes;
    sc++;
  }

  return EXIT_SUCCESS;
}

template <typename Tin, typename T, typename Tout>
int sutra_rtc<Tin, T, Tout>::set_centroids_ref(float *centroids_ref) {
  typename vector<sutra_centroider<Tin, T> *>::iterator sc;
  sc = this->d_centro.begin();
  int inds;
  inds = 0;
  while (sc != this->d_centro.end()) {
    (*sc)->set_centroids_ref(&centroids_ref[inds]);
    inds += (*sc)->nslopes;
    sc++;
  }

  return EXIT_SUCCESS;
}

template <typename Tin, typename T, typename Tout>
int sutra_rtc<Tin, T, Tout>::do_control(int ncntrl) {
  if (this->d_control[ncntrl]->open_loop) {
    return EXIT_SUCCESS;
  } else
    this->d_control[ncntrl]->comp_com();

  return EXIT_SUCCESS;
}

template <typename Tin, typename T, typename Tout>
int sutra_rtc<Tin, T, Tout>::do_clipping(int ncntrl) {
  this->d_control[ncntrl]->clip_commands();

  return EXIT_SUCCESS;
}

template <typename Tin, typename T, typename Tout>
int sutra_rtc<Tin, T, Tout>::comp_voltage(int ncntrl) {
  this->d_control[ncntrl]->comp_voltage();
  // this->d_control[ncntrl]->command_delay();

  return EXIT_SUCCESS;
}

template <typename Tin, typename T, typename Tout>
int sutra_rtc<Tin, T, Tout>::apply_control(int ncntrl, bool compVoltage) {
  if (compVoltage) comp_voltage(ncntrl);

  vector<sutra_dm *>::iterator p;
  p = this->d_control[ncntrl]->d_dmseen.begin();
  int idx = 0;
  // if ((this->d_control[ncntrl]->get_type().compare("ls") == 0) ||
  //     (this->d_control[ncntrl]->get_type().compare("mv") == 0) ||
  //     (this->d_control[ncntrl]->get_type().compare("geo") == 0)) {
  //   // "streamed" controllers case

  //   while (p != this->d_control[ncntrl]->d_dmseen.end()) {
  //     sutra_dm *dm = *p;

  //     int nstreams = this->d_control[ncntrl]->streams->get_nbStreams();
  //     if (nstreams > dm->nactus) {
  //       for (int i = 0; i < nstreams; i++) {
  //         int istart = i * dm->nactus / nstreams;
  //         carmaSafeCall(cudaMemcpyAsync(
  //             dm->d_com->getDataAt(istart),
  //             this->d_control[ncntrl]->d_voltage->getDataAt(idx + istart),
  //             sizeof(float) * dm->nactus / nstreams,
  //             cudaMemcpyDeviceToDevice,
  //             (*this->d_control[ncntrl]->streams)[i]));
  //         dm->comp_shape();
  //       }
  //     } else {
  //       dm->comp_shape(this->d_control[ncntrl]->d_voltage->getDataAt(idx));
  //     }
  //     idx += dm->nactus;
  //     p++;
  //   }
  // } else {  // "non-streamed" controllers
  while (p != this->d_control[ncntrl]->d_dmseen.end()) {
    sutra_dm *dm = *p;
    dm->comp_shape(this->d_control[ncntrl]->d_voltage->getDataAt(idx));
    idx += dm->nactus;
    p++;
  }

  return EXIT_SUCCESS;
}

template class sutra_rtc<float, float, float>;
template class sutra_rtc<float, float, uint16_t>;
template class sutra_rtc<uint16_t, float, float>;
template class sutra_rtc<uint16_t, float, uint16_t>;
#ifdef CAN_DO_HALF
template class sutra_rtc<float, half, float>;
template class sutra_rtc<float, half, uint16_t>;
template class sutra_rtc<uint16_t, half, float>;
template class sutra_rtc<uint16_t, half, uint16_t>;
#endif
