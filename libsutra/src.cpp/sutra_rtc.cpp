// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2019 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.
//  Distributed under GNU - LGPL
//
//  COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
//  General Public License as published by the Free Software Foundation, either version 3 of the License,
//  or any later version.
//
//  COMPASS: End-to-end AO simulation tool using GPU acceleration
//  The COMPASS platform was designed to meet the need of high-performance for the simulation of AO systems.
//
//  The final product includes a software package for simulating all the critical subcomponents of AO,
//  particularly in the context of the ELT and a real-time core based on several control approaches,
//  with performances consistent with its integration into an instrument. Taking advantage of the specific
//  hardware architecture of the GPU, the COMPASS tool allows to achieve adequate execution speeds to
//  conduct large simulation campaigns called to the ELT.
//
//  The COMPASS platform can be used to carry a wide variety of simulations to both testspecific components
//  of AO of the E-ELT (such as wavefront analysis device with a pyramid or elongated Laser star), and
//  various systems configurations such as multi-conjugate AO.
//
//  COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
//  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//  See the GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License along with COMPASS.
//  If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>.
// -----------------------------------------------------------------------------

//! \file      sutra_rtc.cpp
//! \ingroup   libsutra
//! \class     SutraRtc
//! \brief     this class provides the rtc features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.0.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#include <sutra_rtc.h>
#include <sutra_wfs_pyr_pyrhr.h>
#include <sutra_wfs_sh.h>

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
int SutraRtc<Tin, T, Tout>::remove_centroider(int ncentro) {
  delete this->d_centro[ncentro];
  this->d_centro.erase(this->d_centro.begin() + ncentro);
  return EXIT_SUCCESS;
}

template <typename Tin, typename T, typename Tout>
int SutraRtc<Tin, T, Tout>::remove_controller(int ncontrol) {
  delete this->d_control[ncontrol];
  this->d_control.erase(this->d_control.begin() + ncontrol);
  return EXIT_SUCCESS;
}

template <typename Tin, typename T, typename Tout>
int SutraRtc<Tin, T, Tout>::add_centroider(CarmaContext *context, long nvalid,
                                            float offset, float scale,
                                            bool filter_TT, long device,
                                            std::string typec) {
  return add_centroider(context, nvalid, offset, scale, filter_TT, device,
                        typec, nullptr);
}

template <typename Tin, typename T, typename Tout>
int SutraRtc<Tin, T, Tout>::add_centroider(CarmaContext *context, long nvalid,
                                            float offset, float scale,
                                            bool filter_TT, long device,
                                            std::string typec, SutraWfs *wfs) {
  return add_centroider_impl(context, nvalid, offset, scale, filter_TT, device,
                             typec, wfs, std::is_same<T, half>());
}

template <typename Tin, typename T, typename Tout>
template <typename Q>
typename std::enable_if<!std::is_same<Q, half>::value, int>::type
SutraRtc<Tin, T, Tout>::add_centroider_impl(CarmaContext *context,
                                             long nvalid, float offset,
                                             float scale, bool filter_TT,
                                             long device, std::string typec,
                                             SutraWfs *wfs, std::false_type) {
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
          context, wfs, nvalid, 4, offset, scale, filter_TT, device));
    } else if (wfs->type == "pyrhr") {
      SutraWfs_PyrHR *pwfs = dynamic_cast<SutraWfs_PyrHR *>(wfs);
      this->d_centro.push_back(new SutraCentroiderMaskedPix<Tin, T>(
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
int SutraRtc<Tin, T, Tout>::add_centroider_impl(CarmaContext *context,
                                                 long nvalid, float offset,
                                                 float scale, bool filter_TT,
                                                 long device, std::string typec,
                                                 SutraWfs *wfs,
                                                 std::true_type) {
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
          context, wfs, nvalid, 4, offset, scale, filter_TT, device));
    } else if (wfs->type == "pyrhr") {
      SutraWfs_PyrHR *pwfs = dynamic_cast<SutraWfs_PyrHR *>(wfs);
      this->d_centro.push_back(new SutraCentroiderMaskedPix<Tin, T>(
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
int SutraRtc<Tin, T, Tout>::add_controller(CarmaContext *context, int nvalid,
                                            int nslope, int nactu, float delay,
                                            long device, std::string typec,
                                            SutraDms *dms, int *idx_dms,
                                            int ndm,  int *idx_centro, int ncentro, int Nphi,
                                            bool wfs_direction, int nstates) {
  return add_controller_impl(context, this->d_control, nvalid, nslope, nactu,
                             delay, device, typec, dms, idx_dms, ndm, idx_centro, ncentro, Nphi,
                             wfs_direction, nstates, std::is_same<T, half>());
}

template <typename Tin, typename T, typename Tout>
template <typename Q>
typename std::enable_if<!std::is_same<Q, half>::value, int>::type
SutraRtc<Tin, T, Tout>::add_controller_impl(
    CarmaContext *context, vector<SutraController<T, Tout> *> &d_control,
    int nvalid, int nslope, int nactu, float delay, long device,
    std::string typec, SutraDms *dms, int *idx_dms, int ndm,  int *idx_centro, int ncentro, int Nphi,
    bool wfs_direction, int nstates, std::false_type) {
  if (typec.compare("ls") == 0) {
    d_control.push_back(new sutra_controller_ls<T, Tout>(
        context, nvalid, nslope, nactu, delay, dms, idx_dms, ndm, idx_centro, ncentro));
  } else if (typec.compare("geo") == 0) {
    d_control.push_back(new sutra_controller_geo<T, Tout>(
        context, nactu, Nphi, delay, dms, idx_dms, ndm, idx_centro, ncentro, wfs_direction));

  } else if (typec.compare("cured") == 0) {
    d_control.push_back(new SutraControllerCured<T, Tout>(
        context, nvalid, nslope, nactu, delay, dms, idx_dms, ndm, idx_centro, ncentro));
  } else if (typec.compare("mv") == 0) {
    d_control.push_back(new sutra_controller_mv<T, Tout>(
        context, nvalid, nslope, nactu, delay, dms, idx_dms, ndm, idx_centro, ncentro));
  } else if (typec.compare("generic") == 0) {
    d_control.push_back(new sutra_controller_generic<T, Tout>(
        context, nvalid, nslope, nactu, delay, dms, idx_dms, ndm, idx_centro, ncentro, nstates));
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
int SutraRtc<Tin, T, Tout>::add_controller_impl(
    CarmaContext *context, vector<SutraController<T, Tout> *> &d_control,
    int nvalid, int nslope, int nactu, float delay, long device,
    std::string typec, SutraDms *dms, int *idx_dms, int ndm, int *idx_centro, int ncentro, int Nphi,
    bool wfs_direction, int nstates, std::true_type) {
  if (typec.compare("generic") == 0) {
    d_control.push_back(new sutra_controller_generic<T, Tout>(
        context, nvalid, nslope, nactu, delay, dms, idx_dms, ndm, idx_centro, ncentro, nstates));

  } else {
    DEBUG_TRACE("Not implemented in half precision yet");
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

template <typename Tin, typename T, typename Tout>
int SutraRtc<Tin, T, Tout>::do_imat(int ncntrl, SutraDms *ydm) {
  return do_imat_impl(ncntrl, ydm, std::is_same<T, float>());
}

template <typename Tin, typename T, typename Tout>
int SutraRtc<Tin, T, Tout>::do_imat_impl(int ncntrl, SutraDms *ydm,
                                          std::false_type) {
  DEBUG_TRACE("Not implemented for this computation type");
  return EXIT_FAILURE;
}

template <typename Tin, typename T, typename Tout>
template <typename Q>
typename std::enable_if<std::is_same<Q, float>::value, int>::type
SutraRtc<Tin, T, Tout>::do_imat_impl(int ncntrl, SutraDms *ydm,
                                      std::true_type) {
  CarmaObj<T> *d_imat = NULL;
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

  vector<SutraDm *>::iterator p;

  p = this->d_control[ncntrl]->d_dmseen.begin();

  while (p != this->d_control[ncntrl]->d_dmseen.end()) {
    SutraDm *dm = *p;
    if (dm->type == "tt") {
      for (size_t idxh = 0; idxh < this->d_control[ncntrl]->centro_idx.size();
           idxh++) {
        int idx_cntr = this->d_control[ncntrl]->centro_idx[idxh];
        if (this->d_centro[idx_cntr]->filter_TT) {
          std::cout << "Measuring TT reference for centro : " << idx_cntr
                    << std::endl;
          this->d_centro[idx_cntr]->filter_TT = false;

          // Tip Push
          dm->comp_oneactu(0, dm->push4imat);

          this->comp_images_imat(ydm);
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
          this->comp_images_imat(ydm);
          this->d_centro[idx_cntr]->get_cog(
              this->d_centro[idx_cntr]->d_intensities->get_data(),
              this->d_centro[idx_cntr]->d_centro_filtered->get_data(), true);
          float alphai = -.5f / dm->push4imat;

          this->d_centro[idx_cntr]->d_ref_Tip->axpy(
              T(alphai), this->d_centro[idx_cntr]->d_centro_filtered, 1, 1);

          dm->reset_shape();

          // Tilt Push
          dm->comp_oneactu(1, dm->push4imat);

          this->comp_images_imat(ydm);
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
          this->comp_images_imat(ydm);
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
  int inds1 = 0;
  int cc2 = 0;

  std::cout << "Doing imat..." << std::endl;
  while (p != this->d_control[ncntrl]->d_dmseen.end()) {
    SutraDm *dm = *p;
    std::string desc = "DM" + carma_utils::to_string(cc2);
    auto progress = carma_utils::ProgressBar(dm->nactus, desc);
    for (int j = 0; j < dm->nactus; ++j) {
      // Push
      dm->comp_oneactu(j, dm->push4imat);

      this->comp_images_imat(ydm);
      do_centroids(ncntrl, true);

      int device = this->d_control[ncntrl]->d_centroids->get_device();
      this->d_control[ncntrl]->d_centroids->scale(0.5f / dm->push4imat, 1);
      this->d_control[ncntrl]->d_centroids->copy_into(
          d_imat->get_data_at(inds1), this->d_control[ncntrl]->nslope());
      dm->reset_shape();
      // Pull
      dm->comp_oneactu(j, -1.0f * dm->push4imat);
      this->comp_images_imat(ydm);
      device = this->d_control[ncntrl]->d_centroids->get_device();
      do_centroids(ncntrl, true);

      float alphai = -.5f / dm->push4imat;
      cublasSaxpy(this->d_control[ncntrl]->current_context->get_cublas_handle(),
                  this->d_control[ncntrl]->d_centroids->get_nb_elements(), &alphai,
                  this->d_control[ncntrl]->d_centroids->get_data(), 1,
                  d_imat->get_data_at(inds1), 1);

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
int SutraRtc<Tin, T, Tout>::do_imat_basis(int ncntrl, SutraDms *ydm,
                                           int nModes, T *m2v, T *pushAmpl) {
  return do_imat_basis_impl(ncntrl, ydm, nModes, m2v, pushAmpl,
                            std::is_same<T, float>());
}

template <typename Tin, typename T, typename Tout>
int SutraRtc<Tin, T, Tout>::do_imat_basis_impl(int ncntrl, SutraDms *ydm,
                                                int nModes, T *m2v, T *pushAmpl,
                                                std::false_type) {
  DEBUG_TRACE("Not implemented for this computation type");
  return EXIT_FAILURE;
}

template <typename Tin, typename T, typename Tout>
template <typename Q>
typename std::enable_if<std::is_same<Q, float>::value, int>::type
SutraRtc<Tin, T, Tout>::do_imat_basis_impl(int ncntrl, SutraDms *ydm,
                                            int nModes, T *m2v, T *pushAmpl,
                                            std::true_type) {
  CarmaObj<T> d_m2v(this->d_control[ncntrl]->current_context,
                     std::vector<long>{2, ydm->nact_total(), nModes}.data(),
                     m2v);
  CarmaObj<T> d_comm(this->d_control[ncntrl]->current_context,
                      std::vector<long>{1, ydm->nact_total()}.data());

  CarmaObj<T> *d_imat = NULL;

  if (this->d_control[ncntrl]->get_type().compare("ls") == 0) {
    sutra_controller_ls<T, Tout> *control =
        dynamic_cast<sutra_controller_ls<T, Tout> *>(this->d_control[ncntrl]);
    long dims[3] = {2, control->d_imat->get_dims()[1], nModes};
    d_imat = new CarmaObj<T>(this->d_control[ncntrl]->current_context, dims);
    delete control->d_imat;
    control->d_imat = d_imat;
  } else if (this->d_control[ncntrl]->get_type().compare("mv") == 0) {
    sutra_controller_mv<T, Tout> *control =
        dynamic_cast<sutra_controller_mv<T, Tout> *>(this->d_control[ncntrl]);
    long dims[3] = {2, control->d_imat->get_dims()[1], nModes};
    d_imat = new CarmaObj<T>(this->d_control[ncntrl]->current_context, dims);
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
    d_comm.copy_from(d_m2v.get_data_at(j * ydm->nact_total()), ydm->nact_total());
    d_comm.scale(pushAmpl[j], 1);
    // Push
    int actuCount = 0;
    vector<SutraDm *>::iterator p = this->d_control[ncntrl]->d_dmseen.begin();
    while (p != this->d_control[ncntrl]->d_dmseen.end()) {
      // Set each dm
      SutraDm *dm = *p;
      dm->comp_shape(d_comm.get_data_at(actuCount));
      actuCount += dm->nactus;
      ++p;
    }
    this->comp_images_imat(ydm);  // Raytrace & compute all WFS
    do_centroids(ncntrl, true);
    int device = this->d_control[ncntrl]->d_centroids->get_device();
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
    this->comp_images_imat(ydm);  // Raytrace & compute all WFS
    do_centroids(ncntrl, true);

    device = this->d_control[ncntrl]->d_centroids->get_device();
    float alphai = -0.5f / pushAmpl[j];
    cublasSaxpy(this->d_control[ncntrl]->current_context->get_cublas_handle(),
                this->d_control[ncntrl]->d_centroids->get_nb_elements(), &alphai,
                this->d_control[ncntrl]->d_centroids->get_data(), 1,
                d_imat->get_data_at(inds1), 1);

    inds1 += this->d_control[ncntrl]->nslope();
    progress.update();
  }
  progress.finish();
  printf("\n");
  return EXIT_SUCCESS;
}

template <typename Tin, typename T, typename Tout>
int SutraRtc<Tin, T, Tout>::comp_images_imat(SutraDms *ydm) {
  for (size_t idx_cntr = 0; idx_cntr < (this->d_centro).size(); idx_cntr++) {
    SutraWfs *wfs = this->d_centro[idx_cntr]->wfs;
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
int SutraRtc<Tin, T, Tout>::do_imat_geom(int ncntrl, SutraDms *ydm,
                                          int type) {
  return do_imat_geom_impl(ncntrl, ydm, type, std::is_same<T, float>());
}

template <typename Tin, typename T, typename Tout>
int SutraRtc<Tin, T, Tout>::do_imat_geom_impl(int ncntrl, SutraDms *ydm,
                                               int type, std::false_type) {
  DEBUG_TRACE("Not implemented for this computation type");
  return EXIT_FAILURE;
}

template <typename Tin, typename T, typename Tout>
template <typename Q>
typename std::enable_if<std::is_same<Q, float>::value, int>::type
SutraRtc<Tin, T, Tout>::do_imat_geom_impl(int ncntrl, SutraDms *ydm, int type,
                                           std::true_type) {
  vector<SutraDm *>::iterator p;
  p = this->d_control[ncntrl]->d_dmseen.begin();
  int inds1, inds2;
  inds1 = 0;
  while (p != this->d_control[ncntrl]->d_dmseen.end()) {
    SutraDm *dm = *p;
    for (int j = 0; j < dm->nactus; j++) {
      dm->comp_oneactu(j, dm->push4imat);  //
      inds2 = 0;
      for (size_t idx_cntr = 0; idx_cntr < (this->d_control).size();
           idx_cntr++) {
        CarmaObj<T> *d_imat;
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
int SutraRtc<Tin, T, Tout>::do_calibrate_img() {
  for (size_t idx_cntr = 0; idx_cntr < (this->d_centro).size(); idx_cntr++) {
    this->do_calibrate_img(idx_cntr);
  }

  return EXIT_SUCCESS;
}

template <typename Tin, typename T, typename Tout>
int SutraRtc<Tin, T, Tout>::do_calibrate_img(int ncntrl) {
  for (size_t idxh = 0; idxh < this->d_control[ncntrl]->centro_idx.size(); idxh++) {
    int idx_cntr = this->d_control[ncntrl]->centro_idx[idxh];
    this->d_centro[idx_cntr]->calibrate_img();
  }

  return EXIT_SUCCESS;
}

template <typename Tin, typename T, typename Tout>
int SutraRtc<Tin, T, Tout>::do_centroids() {
  for (size_t idx_cntr = 0; idx_cntr < (this->d_centro).size(); idx_cntr++) {
    this->d_centro[idx_cntr]->get_cog();
  }

  return EXIT_SUCCESS;
}

template <typename Tin, typename T, typename Tout>
int SutraRtc<Tin, T, Tout>::do_centroids(int ncntrl) {
  return do_centroids(ncntrl, true);
}

template <typename Tin, typename T, typename Tout>
int SutraRtc<Tin, T, Tout>::do_centroids(int ncntrl, bool noise) {
  int indslope = 0;

  for (size_t idxh = 0; idxh < this->d_control[ncntrl]->centro_idx.size(); idxh++) {
    int idx_cntr = this->d_control[ncntrl]->centro_idx[idxh];
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
int SutraRtc<Tin, T, Tout>::do_centroids_geom(int ncntrl) {
  return do_centroids_geom_impl(ncntrl, std::is_same<T, float>());
}
template <typename Tin, typename T, typename Tout>
int SutraRtc<Tin, T, Tout>::do_centroids_geom_impl(int ncntrl,
                                                    std::false_type) {
  DEBUG_TRACE("Not implemented for this compuation type");
  return EXIT_FAILURE;
}

template <typename Tin, typename T, typename Tout>
template <typename Q>
typename std::enable_if<std::is_same<Q, float>::value, int>::type
SutraRtc<Tin, T, Tout>::do_centroids_geom_impl(int ncntrl, std::true_type) {
  int inds2 = 0;

  for (size_t idxh = 0; idxh < this->d_control[ncntrl]->centro_idx.size(); idxh++) {
    int idx_cntr = this->d_control[ncntrl]->centro_idx[idxh];
    SutraWfs *wfs = this->d_centro[idx_cntr]->wfs;
    if (wfs->type == "sh") {
      SutraWfsSH *_wfs = dynamic_cast<SutraWfsSH *>(wfs);
      _wfs->slopes_geom(this->d_control[ncntrl]->d_centroids->get_data_at(inds2));
    } else if (wfs->type == "pyrhr") {
      SutraWfs_PyrHR *_wfs = dynamic_cast<SutraWfs_PyrHR *>(wfs);
      _wfs->slopes_geom(this->d_control[ncntrl]->d_centroids->get_data_at(inds2));
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
int SutraRtc<Tin, T, Tout>::do_centroids_ref(int ncntrl) {
  typename vector<SutraCentroider<Tin, T> *>::iterator sc;
  sc = this->d_centro.begin();
  while (sc != this->d_centro.end()) {
    (*sc)->d_centroids_ref->reset();
    sc++;
  }
  this->do_centroids(ncntrl);
  sc = this->d_centro.begin();
  int inds;
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
int SutraRtc<Tin, T, Tout>::set_centroids_ref(float *centroids_ref) {
  typename vector<SutraCentroider<Tin, T> *>::iterator sc;
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
int SutraRtc<Tin, T, Tout>::do_control(int ncntrl) {
  if (this->d_control[ncntrl]->open_loop) {
    return EXIT_SUCCESS;
  } else
    this->d_control[ncntrl]->comp_com();

  return EXIT_SUCCESS;
}

template <typename Tin, typename T, typename Tout>
int SutraRtc<Tin, T, Tout>::do_clipping(int ncntrl) {
  this->d_control[ncntrl]->clip_commands();

  return EXIT_SUCCESS;
}

template <typename Tin, typename T, typename Tout>
int SutraRtc<Tin, T, Tout>::comp_voltage(int ncntrl) {
  this->d_control[ncntrl]->comp_voltage();
  // this->d_control[ncntrl]->command_delay();

  return EXIT_SUCCESS;
}

template <typename Tin, typename T, typename Tout>
int SutraRtc<Tin, T, Tout>::apply_control(int ncntrl, bool compVoltage) {
  if (compVoltage) comp_voltage(ncntrl);

  vector<SutraDm *>::iterator p;
  p = this->d_control[ncntrl]->d_dmseen.begin();
  int idx = 0;
  // if ((this->d_control[ncntrl]->get_type().compare("ls") == 0) ||
  //     (this->d_control[ncntrl]->get_type().compare("mv") == 0) ||
  //     (this->d_control[ncntrl]->get_type().compare("geo") == 0)) {
  //   // "streamed" controllers case

  //   while (p != this->d_control[ncntrl]->d_dmseen.end()) {
  //     SutraDm *dm = *p;

  //     int nstreams = this->d_control[ncntrl]->streams->get_nb_streams();
  //     if (nstreams > dm->nactus) {
  //       for (int i = 0; i < nstreams; i++) {
  //         int istart = i * dm->nactus / nstreams;
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
