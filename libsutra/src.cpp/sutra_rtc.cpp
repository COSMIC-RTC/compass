// This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
//
// COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
// General Public License as published by the Free Software Foundation, either version 3 of the 
// License, or any later version.
//
// COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
// See the GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License along with COMPASS. 
// If not, see <https://www.gnu.org/licenses/>
//
//  Copyright (C) 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>

//! \file      sutra_rtc.cpp
//! \ingroup   libsutra
//! \class     SutraRtc
//! \brief     this class provides the rtc features to COMPASS
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \version   5.5.0
//! \date      2022/01/24

#include <sutra_rtc.hpp>
#include <sutra_wfs_pyr_pyrhr.hpp>
#include <sutra_wfs_sh.hpp>

using namespace indicators;

template <typename Tin, typename T, typename Tout>
SutraRtc<Tin, T, Tout>::SutraRtc() {}

template <typename Tin, typename T, typename Tout>
SutraRtc<Tin, T, Tout>::~SutraRtc() {
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
int32_t SutraRtc<Tin, T, Tout>::remove_centroider(int32_t ncentro) {
  delete this->d_centro[ncentro];
  this->d_centro.erase(this->d_centro.begin() + ncentro);
  return EXIT_SUCCESS;
}

template <typename Tin, typename T, typename Tout>
int32_t SutraRtc<Tin, T, Tout>::remove_controller(int32_t ncontrol) {
  delete this->d_control[ncontrol];
  this->d_control.erase(this->d_control.begin() + ncontrol);
  return EXIT_SUCCESS;
}

template <typename Tin, typename T, typename Tout>
int32_t SutraRtc<Tin, T, Tout>::add_centroider(CarmaContext *context, int64_t nvalid,
                                           float offset, float scale,
                                           bool filter_TT, int64_t device,
                                           std::string typec) {
  return add_centroider(context, nvalid, offset, scale, filter_TT, device,
                        typec, nullptr);
}

template <typename Tin, typename T, typename Tout>
int32_t SutraRtc<Tin, T, Tout>::add_centroider(CarmaContext *context, int64_t nvalid,
                                           float offset, float scale,
                                           bool filter_TT, int64_t device,
                                           std::string typec, SutraWfs *wfs) {
  return add_centroider_impl(context, nvalid, offset, scale, filter_TT, device,
                             typec, wfs, std::is_same<T, half>());
}

template <typename Tin, typename T, typename Tout>
template <typename Q>
typename std::enable_if<!std::is_same<Q, half>::value, int32_t>::type
SutraRtc<Tin, T, Tout>::add_centroider_impl(CarmaContext *context, int64_t nvalid,
                                            float offset, float scale,
                                            bool filter_TT, int64_t device,
                                            std::string typec, SutraWfs *wfs,
                                            std::false_type) {
  if (typec.compare("bpcog") == 0)
    this->d_centro.push_back(new SutraCentroiderBpcog<Tin, T>(
        context, wfs, nvalid, offset, scale, filter_TT, device, 10));
  else if (typec.compare("cog") == 0)
    this->d_centro.push_back(new SutraCentroiderCog<Tin, T>(
        context, wfs, nvalid, offset, scale, filter_TT, device));
  else if (typec.compare("corr") == 0)
    this->d_centro.push_back(new SutraCentroiderCorr<Tin, T>(
        context, wfs, nvalid, offset, scale, filter_TT, device));
  else if (typec.compare("pyr") == 0)
    this->d_centro.push_back(new SutraCentroiderPyr<Tin, T>(
        context, wfs, nvalid, offset, scale, filter_TT, device));
  else if (typec.compare("tcog") == 0)
    this->d_centro.push_back(new SutraCentroiderTcog<Tin, T>(
        context, wfs, nvalid, offset, scale, filter_TT, device));
  else if (typec.compare("wcog") == 0)
    this->d_centro.push_back(new SutraCentroiderWcog<Tin, T>(
        context, wfs, nvalid, offset, scale, filter_TT, device));
  else if (typec.compare("maskedpix") == 0) {
    if (wfs == nullptr) {
      this->d_centro.push_back(new SutraCentroiderMaskedPix<Tin, T>(
          context, wfs, nvalid, offset, scale, filter_TT, device));
    } else if (wfs->type == "pyrhr") {
      SutraWfs_PyrHR *pwfs = dynamic_cast<SutraWfs_PyrHR *>(wfs);
      this->d_centro.push_back(new SutraCentroiderMaskedPix<Tin, T>(
          context, pwfs, nvalid * pwfs->npupils, offset, scale, filter_TT,
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
int32_t SutraRtc<Tin, T, Tout>::add_centroider_impl(CarmaContext *context,
                                                int64_t nvalid, float offset,
                                                float scale, bool filter_TT,
                                                int64_t device, std::string typec,
                                                SutraWfs *wfs, std::true_type) {
  if (typec.compare("cog") == 0)
    this->d_centro.push_back(new SutraCentroiderCog<Tin, T>(
        context, wfs, nvalid, offset, scale, filter_TT, device));
  else if (typec.compare("bpcog") == 0)
    this->d_centro.push_back(new SutraCentroiderBpcog<Tin, T>(
        context, wfs, nvalid, offset, scale, filter_TT, device, 10));
  else if (typec.compare("tcog") == 0)
    this->d_centro.push_back(new SutraCentroiderTcog<Tin, T>(
        context, wfs, nvalid, offset, scale, filter_TT, device));
  else if (typec.compare("maskedpix") == 0) {
    if (wfs == nullptr) {
      this->d_centro.push_back(new SutraCentroiderMaskedPix<Tin, T>(
          context, wfs, nvalid, offset, scale, filter_TT, device));
    } else if (wfs->type == "pyrhr") {
      SutraWfs_PyrHR *pwfs = dynamic_cast<SutraWfs_PyrHR *>(wfs);
      this->d_centro.push_back(new SutraCentroiderMaskedPix<Tin, T>(
          context, pwfs, nvalid * pwfs->npupils, offset, scale, filter_TT,
          device));
    } else
      DEBUG_TRACE("WFS must be pyrhr");
  } else
    DEBUG_TRACE("Not implemented for half precision yet");

  return EXIT_SUCCESS;
}
#endif

template <typename Tin, typename T, typename Tout>
int32_t SutraRtc<Tin, T, Tout>::add_controller(CarmaContext *context, std::string typec,int64_t device,
                    float delay, int32_t nslope, int32_t nactu,
                    int32_t nslope_buffers, int32_t nstates, int32_t nstate_buffers, int32_t nmodes,
                    int32_t niir_in, int32_t niir_out, bool polc,bool is_modal,
                    SutraDms *dms, int32_t *idx_dms, int32_t ndm, int32_t *idx_centro, int32_t ncentro,
                    int32_t Nphi, bool wfs_direction) {
  return add_controller_impl(context, this->d_control, typec,device,
                        delay, nslope, nactu, nslope_buffers, nstates, nstate_buffers, nmodes,
                        niir_in, niir_out, polc, is_modal, dms, idx_dms, ndm,
                        idx_centro, ncentro, Nphi, wfs_direction, std::is_same<T, half>());
}

template <typename Tin, typename T, typename Tout>
template <typename Q>
typename std::enable_if<!std::is_same<Q, half>::value, int32_t>::type
SutraRtc<Tin, T, Tout>::add_controller_impl(
    CarmaContext *context, vector<SutraController<T, Tout> *> &d_control,
    std::string typec,int64_t device, float delay, int32_t nslope, int32_t nactu,
    int32_t nslope_buffers, int32_t nstates, int32_t nstate_buffers, int32_t nmodes,
    int32_t niir_in, int32_t niir_out, bool polc,bool is_modal,
    SutraDms *dms, int32_t *idx_dms, int32_t ndm, int32_t *idx_centro, int32_t ncentro,
    int32_t Nphi, bool wfs_direction, std::false_type) {
  if (typec.compare("ls") == 0) {
    d_control.push_back(new SutraControllerLs<T, Tout>(
        context, nslope, nactu, delay, dms, idx_dms, ndm, idx_centro,
        ncentro));
  } else if (typec.compare("geo") == 0) {
    d_control.push_back(new SutraControllerGeo<T, Tout>(
        context, nactu, Nphi, delay, dms, idx_dms, ndm, idx_centro, ncentro,
        wfs_direction));

  } else if (typec.compare("mv") == 0) {
    d_control.push_back(new SutraControllerMv<T, Tout>(
        context, nslope, nactu, delay, dms, idx_dms, ndm, idx_centro,
        ncentro));
  } else if (typec.compare("generic") == 0) {
    d_control.push_back(new SutraControllerGeneric<T, Tout>(
        context, nslope, nactu, delay, dms, idx_dms, ndm, idx_centro,
        ncentro, nstates));
  } else if(typec.compare("generic_linear") == 0){
    d_control.push_back(new SutraControllerGenericLinear<T, Tout>( context,
    nslope, nslope_buffers, nactu, nstates, nstate_buffers, nmodes,
    niir_in, niir_out, delay, polc, is_modal,
    dms, idx_dms, ndm, idx_centro, ncentro));
  } else {
    DEBUG_TRACE("Controller '%s' unknown\n", typec.c_str());
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

template <typename Tin, typename T, typename Tout>
int32_t SutraRtc<Tin, T, Tout>::add_controller_impl(
    CarmaContext *context, vector<SutraController<T, Tout> *> &d_control,
    std::string typec,int64_t device, float delay, int32_t nslope, int32_t nactu,
    int32_t nslope_buffers, int32_t nstates, int32_t nstate_buffers, int32_t nmodes,
    int32_t niir_in, int32_t niir_out, bool polc,bool is_modal,
    SutraDms *dms, int32_t *idx_dms, int32_t ndm, int32_t *idx_centro, int32_t ncentro,
    int32_t Nphi, bool wfs_direction, std::true_type) {
  if (typec.compare("generic") == 0) {
    d_control.push_back(new SutraControllerGeneric<T, Tout>(
        context, nslope, nactu, delay, dms, idx_dms, ndm, idx_centro,
        ncentro, nstates));

  } else {
    DEBUG_TRACE("Not implemented in half precision yet");
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

template <typename Tin, typename T, typename Tout>
int32_t SutraRtc<Tin, T, Tout>::do_imat(int32_t ncntrl, SutraDms *ydm, int32_t kernconv) {
  return do_imat_impl(ncntrl, ydm, kernconv, std::is_same<T, float>());
}

template <typename Tin, typename T, typename Tout>
int32_t SutraRtc<Tin, T, Tout>::do_imat_impl(int32_t ncntrl, SutraDms *ydm, int32_t kernconv,
                                         std::false_type) {
  DEBUG_TRACE("Not implemented for this computation type");
  return EXIT_FAILURE;
}

template <typename Tin, typename T, typename Tout>
template <typename Q>
typename std::enable_if<std::is_same<Q, float>::value, int32_t>::type
SutraRtc<Tin, T, Tout>::do_imat_impl(int32_t ncntrl, SutraDms *ydm, int32_t kernconv,
                                     std::true_type) {
  CarmaObj<T> *d_imat = NULL;
  if (this->d_control[ncntrl]->get_type().compare("ls") == 0) {
    SutraControllerLs<T, Tout> *control =
        dynamic_cast<SutraControllerLs<T, Tout> *>(this->d_control[ncntrl]);
    d_imat = control->d_imat;
  } else if (this->d_control[ncntrl]->get_type().compare("mv") == 0) {
    SutraControllerMv<T, Tout> *control =
        dynamic_cast<SutraControllerMv<T, Tout> *>(this->d_control[ncntrl]);
    d_imat = control->d_imat;
  } else {
    DEBUG_TRACE("Controller needs to be LS or MV\n");
    return EXIT_SUCCESS;
  }

  vector<SutraDm *>::iterator p;

  p = this->d_control[ncntrl]->d_dmseen.begin();

  while (p != this->d_control[ncntrl]->d_dmseen.end()) {
    SutraDm *dm = *p;
    if (dm->type == "tt") {
      for (size_t idxh = 0; idxh < this->d_control[ncntrl]->centro_idx.size();
           idxh++) {
        int32_t idx_cntr = this->d_control[ncntrl]->centro_idx[idxh];
        if (this->d_centro[idx_cntr]->filter_TT) {
          std::cout << "Measuring TT reference for centro : " << idx_cntr
                    << std::endl;
          this->d_centro[idx_cntr]->filter_TT = false;

          // Tip Push
          dm->comp_oneactu(0, dm->push4imat);

          this->comp_images_imat(ydm, kernconv);
          this->d_centro[idx_cntr]->get_cog(
              this->d_centro[idx_cntr]->d_intensities->get_data(),
              this->d_centro[idx_cntr]->d_centro_filtered->get_data(), true);
          this->d_centro[idx_cntr]->d_centro_filtered->scale(
              0.5f / dm->push4imat, 1);
          this->d_centro[idx_cntr]->d_centro_filtered->copy_into(
              this->d_centro[idx_cntr]->d_ref_Tip->get_data(),
              this->d_centro[idx_cntr]->nslopes);

          dm->reset_shape();

          // Tip Pull
          dm->comp_oneactu(0, -1.0f * dm->push4imat);
          this->comp_images_imat(ydm, kernconv);
          this->d_centro[idx_cntr]->get_cog(
              this->d_centro[idx_cntr]->d_intensities->get_data(),
              this->d_centro[idx_cntr]->d_centro_filtered->get_data(), true);
          float alphai = -.5f / dm->push4imat;

          this->d_centro[idx_cntr]->d_ref_Tip->axpy(
              T(alphai), this->d_centro[idx_cntr]->d_centro_filtered, 1, 1);

          dm->reset_shape();

          // Tilt Push
          dm->comp_oneactu(1, dm->push4imat);

          this->comp_images_imat(ydm, kernconv);
          this->d_centro[idx_cntr]->get_cog(
              this->d_centro[idx_cntr]->d_intensities->get_data(),
              this->d_centro[idx_cntr]->d_centro_filtered->get_data(), true);
          this->d_centro[idx_cntr]->d_centro_filtered->scale(
              0.5f / dm->push4imat, 1);
          this->d_centro[idx_cntr]->d_centro_filtered->copy_into(
              this->d_centro[idx_cntr]->d_ref_Tilt->get_data(),
              this->d_centro[idx_cntr]->nslopes);

          dm->reset_shape();

          // Tilt Pull
          dm->comp_oneactu(1, -1.0f * dm->push4imat);
          this->comp_images_imat(ydm, kernconv);
          this->d_centro[idx_cntr]->get_cog(
              this->d_centro[idx_cntr]->d_intensities->get_data(),
              this->d_centro[idx_cntr]->d_centro_filtered->get_data(), true);

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
  int32_t inds1 = 0;
  int32_t cc2 = 0;

  std::cout << "Doing imat..." << std::endl;
  while (p != this->d_control[ncntrl]->d_dmseen.end()) {
    SutraDm *dm = *p;
    ProgressBar bar{option::BarWidth{50}, option::ForegroundColor{Color::white},
                    option::ShowElapsedTime{true}, option::ShowRemainingTime{true},
                    option::PrefixText{"DM" + carma_utils::to_string(cc2)}, option::MaxProgress{dm->nactus}};
    for (int32_t j = 0; j < dm->nactus; ++j) {
      // Push
      dm->comp_oneactu(j, dm->push4imat);

      this->comp_images_imat(ydm, kernconv);
      do_centroids(ncntrl, true);

      int32_t device = this->d_control[ncntrl]->d_centroids->get_device();
      this->d_control[ncntrl]->d_centroids->scale(0.5f / dm->push4imat, 1);
      this->d_control[ncntrl]->d_centroids->copy_into(
          d_imat->get_data_at(inds1), this->d_control[ncntrl]->nslope());
      dm->reset_shape();
      // Pull
      dm->comp_oneactu(j, -1.0f * dm->push4imat);
      this->comp_images_imat(ydm, kernconv);
      device = this->d_control[ncntrl]->d_centroids->get_device();
      do_centroids(ncntrl, true);

      float alphai = -.5f / dm->push4imat;
      cublasSaxpy(this->d_control[ncntrl]->current_context->get_cublas_handle(),
                  this->d_control[ncntrl]->d_centroids->get_nb_elements(),
                  &alphai, this->d_control[ncntrl]->d_centroids->get_data(), 1,
                  d_imat->get_data_at(inds1), 1);

      dm->reset_shape();
      inds1 += this->d_control[ncntrl]->nslope();
      bar.tick();
      // printf("\rDoing imat...%d%%",(cc*100/nactu));
    }
    ++p;
    ++cc2;
  }
  printf("\n");
  return EXIT_SUCCESS;
}

template <typename Tin, typename T, typename Tout>
int32_t SutraRtc<Tin, T, Tout>::do_imat_basis(int32_t ncntrl, SutraDms *ydm, int32_t nModes,
                                          T *m2v, T *pushAmpl, int32_t kernconv) {
  return do_imat_basis_impl(ncntrl, ydm, nModes, m2v, pushAmpl, kernconv,
                            std::is_same<T, float>());
}

template <typename Tin, typename T, typename Tout>
int32_t SutraRtc<Tin, T, Tout>::do_imat_basis_impl(int32_t ncntrl, SutraDms *ydm,
                                               int32_t nModes, T *m2v, T *pushAmpl, int32_t kernconv,
                                               std::false_type) {
  DEBUG_TRACE("Not implemented for this computation type");
  return EXIT_FAILURE;
}

template <typename Tin, typename T, typename Tout>
template <typename Q>
typename std::enable_if<std::is_same<Q, float>::value, int32_t>::type
SutraRtc<Tin, T, Tout>::do_imat_basis_impl(int32_t ncntrl, SutraDms *ydm,
                                           int32_t nModes, T *m2v, T *pushAmpl, int32_t kernconv,
                                           std::true_type) {
  CarmaObj<T> d_m2v(this->d_control[ncntrl]->current_context,
                    std::vector<int64_t>{2, ydm->nact_total(), nModes}.data(),
                    m2v);
  CarmaObj<T> d_comm(this->d_control[ncntrl]->current_context,
                     std::vector<int64_t>{1, ydm->nact_total()}.data());

  CarmaObj<T> *d_imat = NULL;

  if (this->d_control[ncntrl]->get_type().compare("ls") == 0) {
    SutraControllerLs<T, Tout> *control =
        dynamic_cast<SutraControllerLs<T, Tout> *>(this->d_control[ncntrl]);
    int64_t dims[3] = {2, control->d_imat->get_dims()[1], nModes};
    d_imat = new CarmaObj<T>(this->d_control[ncntrl]->current_context, dims);
    delete control->d_imat;
    control->d_imat = d_imat;
  } else if (this->d_control[ncntrl]->get_type().compare("mv") == 0) {
    SutraControllerMv<T, Tout> *control =
        dynamic_cast<SutraControllerMv<T, Tout> *>(this->d_control[ncntrl]);
    int64_t dims[3] = {2, control->d_imat->get_dims()[1], nModes};
    d_imat = new CarmaObj<T>(this->d_control[ncntrl]->current_context, dims);
    delete control->d_imat;
    control->d_imat = d_imat;
  } else {
    DEBUG_TRACE("Controller needs to be LS or MV\n");
    return EXIT_SUCCESS;
  }
  int32_t inds1 = 0;

  std::cout << "Doing imat modal..." << std::endl;
    ProgressBar bar{option::BarWidth{50}, option::ForegroundColor{Color::white},
                    option::ShowElapsedTime{true}, option::ShowRemainingTime{true},
                    option::PrefixText{"Modal iMat"}, option::MaxProgress{nModes}};
  for (int32_t j = 0; j < nModes; ++j) {
    // For each mode
    d_comm.copy_from(d_m2v.get_data_at(j * ydm->nact_total()),
                     ydm->nact_total());
    d_comm.scale(pushAmpl[j], 1);
    // Push
    int32_t actuCount = 0;
    vector<SutraDm *>::iterator p = this->d_control[ncntrl]->d_dmseen.begin();
    while (p != this->d_control[ncntrl]->d_dmseen.end()) {
      // Set each dm
      SutraDm *dm = *p;
      dm->comp_shape(d_comm.get_data_at(actuCount));
      actuCount += dm->nactus;
      ++p;
    }
    this->comp_images_imat(ydm, kernconv);  // Raytrace & compute all WFS
    do_centroids(ncntrl, true);
    int32_t device = this->d_control[ncntrl]->d_centroids->get_device();
    this->d_control[ncntrl]->d_centroids->scale(0.5f / pushAmpl[j], 1);
    this->d_control[ncntrl]->d_centroids->copy_into(
        d_imat->get_data_at(inds1), this->d_control[ncntrl]->nslope());

    // Pull
    actuCount = 0;
    d_comm.scale(-1.0f, 1);
    p = this->d_control[ncntrl]->d_dmseen.begin();
    while (p != this->d_control[ncntrl]->d_dmseen.end()) {
      // Set each dm
      SutraDm *dm = *p;
      dm->comp_shape(d_comm.get_data_at(actuCount));
      actuCount += dm->nactus;
      ++p;
    }
    this->comp_images_imat(ydm, kernconv);  // Raytrace & compute all WFS
    do_centroids(ncntrl, true);

    device = this->d_control[ncntrl]->d_centroids->get_device();
    float alphai = -0.5f / pushAmpl[j];
    cublasSaxpy(this->d_control[ncntrl]->current_context->get_cublas_handle(),
                this->d_control[ncntrl]->d_centroids->get_nb_elements(),
                &alphai, this->d_control[ncntrl]->d_centroids->get_data(), 1,
                d_imat->get_data_at(inds1), 1);

    inds1 += this->d_control[ncntrl]->nslope();
    bar.tick();
  }
  printf("\n");
  return EXIT_SUCCESS;
}

template <typename Tin, typename T, typename Tout>
int32_t SutraRtc<Tin, T, Tout>::comp_images_imat(SutraDms *ydm, int32_t kernconv) {
  for (size_t idx_cntr = 0; idx_cntr < (this->d_centro).size(); idx_cntr++) {
    SutraWfs *wfs = this->d_centro[idx_cntr]->wfs;
    float tmp_noise = wfs->noise;
    float tmp_nphot = wfs->nphot;
    wfs->nphot = wfs->nphot4imat;
    wfs->noise = -1;
    wfs->kernconv = kernconv;
    wfs->sensor_trace(ydm, 1);
    wfs->comp_image();
    wfs->noise = tmp_noise;
    wfs->nphot = tmp_nphot;
    wfs->kernconv = false;
  }
  return EXIT_SUCCESS;
}

template <typename Tin, typename T, typename Tout>
int32_t SutraRtc<Tin, T, Tout>::do_imat_geom(int32_t ncntrl, SutraDms *ydm, int32_t type) {
  return do_imat_geom_impl(ncntrl, ydm, type, std::is_same<T, float>());
}

template <typename Tin, typename T, typename Tout>
int32_t SutraRtc<Tin, T, Tout>::do_imat_geom_impl(int32_t ncntrl, SutraDms *ydm,
                                              int32_t type, std::false_type) {
  DEBUG_TRACE("Not implemented for this computation type");
  return EXIT_FAILURE;
}

template <typename Tin, typename T, typename Tout>
template <typename Q>
typename std::enable_if<std::is_same<Q, float>::value, int32_t>::type
SutraRtc<Tin, T, Tout>::do_imat_geom_impl(int32_t ncntrl, SutraDms *ydm, int32_t type,
                                          std::true_type) {
  vector<SutraDm *>::iterator p;
  p = this->d_control[ncntrl]->d_dmseen.begin();
  int32_t inds1, inds2;
  inds1 = 0;
  while (p != this->d_control[ncntrl]->d_dmseen.end()) {
    SutraDm *dm = *p;
    for (int32_t j = 0; j < dm->nactus; j++) {
      dm->comp_oneactu(j, dm->push4imat);  //
      inds2 = 0;
      for (size_t idx_cntr = 0; idx_cntr < (this->d_control).size();
           idx_cntr++) {
        CarmaObj<T> *d_imat;
        if (this->d_control[ncntrl]->get_type().compare("ls") == 0) {
          SutraControllerLs<T, Tout> *control =
              dynamic_cast<SutraControllerLs<T, Tout> *>(
                  this->d_control[ncntrl]);
          d_imat = control->d_imat;
        } else if (this->d_control[ncntrl]->get_type().compare("mv") == 0) {
          SutraControllerMv<T, Tout> *control =
              dynamic_cast<SutraControllerMv<T, Tout> *>(
                  this->d_control[ncntrl]);
          d_imat = control->d_imat;
        } else {
          DEBUG_TRACE("controller should be ls or mv");
          return EXIT_FAILURE;
        }

        SutraWfs *wfs = this->d_centro[idx_cntr]->wfs;

        wfs->sensor_trace(ydm, 1);
        // sensors->d_wfs[nwfs]->comp_image();
        if (wfs->type != "sh") {
          SutraWfsSH *_wfs = dynamic_cast<SutraWfsSH *>(wfs);
          _wfs->slopes_geom(d_imat[inds1 + inds2], type);
        } else if (wfs->type != "pyrhr") {
          SutraWfs_PyrHR *_wfs = dynamic_cast<SutraWfs_PyrHR *>(wfs);
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
int32_t SutraRtc<Tin, T, Tout>::do_calibrate_img() {
  for (size_t idx_cntr = 0; idx_cntr < (this->d_centro).size(); idx_cntr++) {
    this->do_calibrate_img(idx_cntr);
  }

  return EXIT_SUCCESS;
}

template <typename Tin, typename T, typename Tout>
int32_t SutraRtc<Tin, T, Tout>::do_calibrate_img(int32_t ncntrl) {
  for (size_t idxh = 0; idxh < this->d_control[ncntrl]->centro_idx.size();
       idxh++) {
    int32_t idx_cntr = this->d_control[ncntrl]->centro_idx[idxh];
    this->d_centro[idx_cntr]->calibrate_img();
  }

  return EXIT_SUCCESS;
}

template <typename Tin, typename T, typename Tout>
int32_t SutraRtc<Tin, T, Tout>::do_centroids() {
  for (size_t idx_cntr = 0; idx_cntr < (this->d_centro).size(); idx_cntr++) {
    this->d_centro[idx_cntr]->get_cog();
  }

  return EXIT_SUCCESS;
}

template <typename Tin, typename T, typename Tout>
int32_t SutraRtc<Tin, T, Tout>::do_centroids(int32_t ncntrl) {
  return do_centroids(ncntrl, true);
}

template <typename Tin, typename T, typename Tout>
int32_t SutraRtc<Tin, T, Tout>::do_centroids(int32_t ncntrl, bool noise) {
  int32_t indslope = 0;

  for (size_t idxh = 0; idxh < this->d_control[ncntrl]->centro_idx.size();
       idxh++) {
    int32_t idx_cntr = this->d_control[ncntrl]->centro_idx[idxh];
    if (this->d_centro[idx_cntr]->wfs != nullptr) {
      this->d_centro[idx_cntr]->get_cog(
          this->d_centro[idx_cntr]->d_intensities->get_data(),
          this->d_control[ncntrl]->d_centroids->get_data_at(indslope), noise);

      indslope += this->d_centro[idx_cntr]->nslopes;
    } else {
      this->d_centro[idx_cntr]->get_cog(
          this->d_centro[idx_cntr]->d_img->get_data(),
          this->d_centro[idx_cntr]->d_intensities->get_data(),
          this->d_control[ncntrl]->d_centroids->get_data_at(indslope),
          this->d_centro[idx_cntr]->nvalid, this->d_centro[idx_cntr]->npix,
          this->d_centro[idx_cntr]->d_img->get_dims(1));

      indslope += this->d_centro[idx_cntr]->nslopes;
    }
  }
  // remove_ref(ncntrl);

  return EXIT_SUCCESS;
}

template <typename Tin, typename T, typename Tout>
int32_t SutraRtc<Tin, T, Tout>::do_centroids_geom(int32_t ncntrl, int32_t type) {
  return do_centroids_geom_impl(ncntrl, type, std::is_same<T, float>());
}

template <typename Tin, typename T, typename Tout>
int32_t SutraRtc<Tin, T, Tout>::do_centroids_geom_impl(int32_t ncntrl, int32_t type,
                                                   std::false_type) {
  DEBUG_TRACE("Not implemented for this compuation type");
  return EXIT_FAILURE;
}

template <typename Tin, typename T, typename Tout>
template <typename Q>
typename std::enable_if<std::is_same<Q, float>::value, int32_t>::type
SutraRtc<Tin, T, Tout>::do_centroids_geom_impl(int32_t ncntrl, int32_t type,
                                               std::true_type) {
  int32_t inds2 = 0;

  for (size_t idxh = 0; idxh < this->d_control[ncntrl]->centro_idx.size();
       idxh++) {
    int32_t idx_cntr = this->d_control[ncntrl]->centro_idx[idxh];
    SutraWfs *wfs = this->d_centro[idx_cntr]->wfs;
    if (wfs->type == "sh") {
      SutraWfsSH *_wfs = dynamic_cast<SutraWfsSH *>(wfs);
      _wfs->slopes_geom(
          this->d_control[ncntrl]->d_centroids->get_data_at(inds2), type);
    } else if (wfs->type == "pyrhr") {
      SutraWfs_PyrHR *_wfs = dynamic_cast<SutraWfs_PyrHR *>(wfs);
      _wfs->slopes_geom(
          this->d_control[ncntrl]->d_centroids->get_data_at(inds2), type);
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
int32_t SutraRtc<Tin, T, Tout>::do_centroids_ref(int32_t ncntrl) {
  typename vector<SutraCentroider<Tin, T> *>::iterator sc;
  sc = this->d_centro.begin();
  while (sc != this->d_centro.end()) {
    (*sc)->d_centroids_ref->reset();
    sc++;
  }
  this->do_centroids(ncntrl);
  sc = this->d_centro.begin();
  int32_t inds;
  inds = 0;
  while (sc != this->d_centro.end()) {
    (*sc)->d_centroids_ref->axpy(1.0f, this->d_control[ncntrl]->d_centroids, 1,
                                 1, inds);
    inds += (*sc)->nslopes;
    sc++;
  }

  return EXIT_SUCCESS;
}

template <typename Tin, typename T, typename Tout>
int32_t SutraRtc<Tin, T, Tout>::set_centroids_ref(float *centroids_ref) {
  typename vector<SutraCentroider<Tin, T> *>::iterator sc;
  sc = this->d_centro.begin();
  int32_t inds;
  inds = 0;
  while (sc != this->d_centro.end()) {
    (*sc)->set_centroids_ref(&centroids_ref[inds]);
    inds += (*sc)->nslopes;
    sc++;
  }

  return EXIT_SUCCESS;
}

template <typename Tin, typename T, typename Tout>
int32_t SutraRtc<Tin, T, Tout>::do_control(int32_t ncntrl) {
  if (this->d_control[ncntrl]->open_loop) {
    return EXIT_SUCCESS;
  } else
    this->d_control[ncntrl]->comp_com();

  return EXIT_SUCCESS;
}

template <typename Tin, typename T, typename Tout>
int32_t SutraRtc<Tin, T, Tout>::do_clipping(int32_t ncntrl) {
  this->d_control[ncntrl]->clip_commands();

  return EXIT_SUCCESS;
}

template <typename Tin, typename T, typename Tout>
int32_t SutraRtc<Tin, T, Tout>::comp_voltage(int32_t ncntrl) {
  this->d_control[ncntrl]->comp_voltage();
  // this->d_control[ncntrl]->command_delay();

  return EXIT_SUCCESS;
}

template <typename Tin, typename T, typename Tout>
int32_t SutraRtc<Tin, T, Tout>::apply_control(int32_t ncntrl, bool compVoltage) {
  if (compVoltage) comp_voltage(ncntrl);

  vector<SutraDm *>::iterator p;
  p = this->d_control[ncntrl]->d_dmseen.begin();
  int32_t idx = 0;
  // if ((this->d_control[ncntrl]->get_type().compare("ls") == 0) ||
  //     (this->d_control[ncntrl]->get_type().compare("mv") == 0) ||
  //     (this->d_control[ncntrl]->get_type().compare("geo") == 0)) {
  //   // "streamed" controllers case

  //   while (p != this->d_control[ncntrl]->d_dmseen.end()) {
  //     SutraDm *dm = *p;

  //     int32_t nstreams = this->d_control[ncntrl]->streams->get_nb_streams();
  //     if (nstreams > dm->nactus) {
  //       for (int32_t i = 0; i < nstreams; i++) {
  //         int32_t istart = i * dm->nactus / nstreams;
  //         carma_safe_call(cudaMemcpyAsync(
  //             dm->d_com->get_data_at(istart),
  //             this->d_control[ncntrl]->d_voltage->get_data_at(idx + istart),
  //             sizeof(float) * dm->nactus / nstreams,
  //             cudaMemcpyDeviceToDevice,
  //             (*this->d_control[ncntrl]->streams)[i]));
  //         dm->comp_shape();
  //       }
  //     } else {
  //       dm->comp_shape(this->d_control[ncntrl]->d_voltage->get_data_at(idx));
  //     }
  //     idx += dm->nactus;
  //     p++;
  //   }
  // } else {  // "non-streamed" controllers
  while (p != this->d_control[ncntrl]->d_dmseen.end()) {
    SutraDm *dm = *p;
    dm->comp_shape(this->d_control[ncntrl]->d_voltage->get_data_at(idx));
    idx += dm->nactus;
    p++;
  }

  return EXIT_SUCCESS;
}

template class SutraRtc<float, float, float>;
template class SutraRtc<float, float, uint16_t>;
template class SutraRtc<uint16_t, float, float>;
template class SutraRtc<uint16_t, float, uint16_t>;
#ifdef CAN_DO_HALF
template class SutraRtc<float, half, float>;
template class SutraRtc<float, half, uint16_t>;
template class SutraRtc<uint16_t, half, float>;
template class SutraRtc<uint16_t, half, uint16_t>;
#endif
