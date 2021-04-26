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

//! \file      sutra_source.cpp
//! \ingroup   libsutra
//! \class     SutraSource
//! \brief     this class provides the source features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.1.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#include <sutra_source.h>

SutraSource::SutraSource(CarmaContext *context, float xpos, float ypos,
                           float lambda, float mag, float zerop, long size,
                           string type, CarmaObj<float> *pupil, int Npts,
                           int device)
    : current_context(context),
      device(device),
      d_pupil(pupil),
      d_ncpa_phase(nullptr) {
  current_context->set_active_device(device, 1);

  this->init_source(context, xpos, ypos, lambda, mag, zerop, size, type,
                    device);
  // float h_pupil[this->d_pupil->get_nb_elements()];
  float *h_pupil = new float[this->d_pupil->get_nb_elements()];
  this->d_pupil->device2host(h_pupil);

  long *dims_data1 = new long[2];
  dims_data1[0] = 1;
  dims_data1[1] = pow(2, (long)(logf(Npts) / logf(2)) + 1);
  // this->d_phasepts = new CarmaObj<float>(this->current_context, dims_data1);
  dims_data1[1] = Npts;
  this->d_phasepts = new CarmaObj<float>(this->current_context, dims_data1);
  this->d_wherephase = new CarmaObj<int>(this->current_context, dims_data1);
  int *wherephase = new int[Npts];
  int cpt = 0;
  for (int cc = 0; cc < this->d_pupil->get_nb_elements(); cc++) {
    if (h_pupil[cc] > 0) {
      wherephase[cpt] = cc;
      cpt += 1;
    }
  }

  this->d_wherephase->host2device(wherephase);
  delete[] wherephase;
  delete[] dims_data1;
  delete[] h_pupil;
}

SutraSource::SutraSource(CarmaContext *context, float xpos, float ypos,
                           float lambda, float mag, float zerop, long size,
                           string type, int device)
    : current_context(context), device(device), d_ncpa_phase(nullptr) {
  this->device = device;
  current_context = context;
  current_context->set_active_device(device, 1);
  this->init_source(context, xpos, ypos, lambda, mag, zerop, size, type,
                    device);
}

inline int SutraSource::init_source(CarmaContext *context, float xpos,
                                     float ypos, float lambda, float mag,
                                     float zerop, long size, string type,
                                     int device) {
  current_context->set_active_device(device, 1);
  this->current_context = context;
  this->strehl_counter = 0;

  this->posx = xpos;
  this->posy = ypos;
  this->lambda = lambda;
  this->mag = mag;
  this->zp = zerop;
  this->npts = size;

  this->d_phase = new SutraPhase(context, size);

  // ADDING INSTRUMENTAL PHASE
  // this->d_phase_instru = new SutraPhase(context, size);
  //

  this->type = type;
  this->device = device;
  this->scale = float(2 * CARMA_PI / lambda);  // phase is expected in microns

  // cudaDeviceProp deviceProperties =
  // current_context->get_device(device)->get_properties();
  // this->blockSize = (int)sqrt(deviceProperties.maxThreadsPerBlock);
  this->block_size = 8;
  this->phase_var_avg = 0;
  this->phase_var_count = -10;

  long *dims_data2 = new long[3];
  dims_data2[0] = 2;
  dims_data2[1] = size;
  dims_data2[2] = size;
  int nstreams =
      (size + this->block_size - size % this->block_size) / this->block_size;

  this->phase_telemetry =
      new CarmaHostObj<float>(dims_data2, MA_WRICOMB, nstreams);

  this->d_image_se = 0L;
  this->d_amplipup = 0L;
  this->d_image_le = 0L;
  this->d_phasepts = 0L;
  this->d_wherephase = 0L;

  if (type != "wfs") {
    int mradix = 2;  // fft_goodsize(size);

    int fft_size = pow(mradix, (long)(logf(2 * size) / logf(mradix)) + 1);
    dims_data2[1] = fft_size;
    dims_data2[2] = fft_size;

    this->d_image_se = new CarmaObj<float>(context, dims_data2);
    this->d_amplipup = new CarmaObj<cuFloatComplex>(context, dims_data2);

    cufftHandle *plan = this->d_amplipup->get_plan();  ///< FFT plan
    carmafft_safe_call(cufftPlan2d(plan, this->d_amplipup->get_dims(1),
                                 this->d_amplipup->get_dims(2), CUFFT_C2C));
  }

  this->lgs = false;
  this->G = 1.0f;
  this->thetaML = 0.0f;
  this->dx = 0.0f;
  this->dy = 0.0f;

  /// temporary array for accurate strehl computation
  dims_data2[1] = 2 * this->d_smallimg_size + 1;
  dims_data2[2] = 2 * this->d_smallimg_size + 1;
  this->d_smallimg = new CarmaObj<float>(context, dims_data2);

  delete[] dims_data2;

  return EXIT_SUCCESS;
}

SutraSource::~SutraSource() {
  // delete this->current_context;
  this->current_context->set_active_device(this->device, 1);
  delete this->d_smallimg;
  delete this->d_phase;
  delete this->phase_telemetry;

  if (this->d_image_se != 0L) delete this->d_image_se;
  if (this->d_image_le != 0L) delete this->d_image_le;
  if (this->d_amplipup != 0L) delete this->d_amplipup;
  if (this->d_phasepts != 0L) delete this->d_phasepts;
  if (this->d_wherephase != 0L) delete this->d_wherephase;
  if (this->d_ncpa_phase != nullptr) delete this->d_ncpa_phase;

  this->xoff.clear();
  this->yoff.clear();
}

int SutraSource::init_strehlmeter() {
  current_context->set_active_device(device, 1);
  this->strehl_counter = 0;
  this->comp_image(1, false);

  cudaMemcpy(&(this->ref_strehl),
             &(this->d_image_se->get_data()[this->d_image_se->aimax(1)]),
             sizeof(float), cudaMemcpyDeviceToHost);

  if (this->d_image_le == 0L)
    this->d_image_le = new CarmaObj<float>(this->current_context,
                                            this->d_image_se->get_dims());
  else {  // Reset strehl case
    carma_safe_call(cudaMemset(this->d_image_le->get_data(), 0,
                             sizeof(float) * this->d_image_le->get_nb_elements()));
    this->strehl_counter = 0;
    this->phase_var_avg = 0.f;
    this->phase_var_count = 0;
  }

  return EXIT_SUCCESS;
}

int SutraSource::reset_strehlmeter() {
  carma_safe_call(cudaMemset(this->d_image_le->get_data(), 0,
                           sizeof(float) * this->d_image_le->get_nb_elements()));
  this->strehl_counter = 0;

  return EXIT_SUCCESS;
}

int SutraSource::reset_phase() {
  this->d_phase->d_screen->reset();
  return EXIT_SUCCESS;
}
int SutraSource::add_layer(string type, int idx, float mxoff, float myoff) {
  current_context->set_active_device(device, 1);
  xoff[make_pair(type, idx)] = mxoff;
  yoff[make_pair(type, idx)] = myoff;

  return EXIT_SUCCESS;
}

int SutraSource::remove_layer(string type, int idx) {
  current_context->set_active_device(device, 1);
  xoff.erase(make_pair(type, idx));
  yoff.erase(make_pair(type, idx));

  return EXIT_SUCCESS;
}

int SutraSource::raytrace(SutraAtmos *yatmos, bool async) {
  //  carma_safe_call(cudaDeviceSynchronize());
  current_context->set_active_device(device, 1);
  carma_safe_call(
      cudaMemset(this->d_phase->d_screen->get_data(), 0,
                 sizeof(float) * this->d_phase->d_screen->get_nb_elements()));

  float delta = 1.0f;;
  map<type_screen, float>::iterator p;
  p = xoff.begin();

  while (p != xoff.end()) {
    string types = p->first.first;

    if (types.find("atmos") == 0) {
      int idx = p->first.second;
      SutraTurbuScreen *ps;
      ps = yatmos->d_screens[idx];
      if ((p == xoff.end()) && async) {
        target_raytrace_async(
            this->phase_telemetry, this->d_phase->d_screen->get_data(),
            ps->d_tscreen->d_screen->get_data(),
            (int)d_phase->d_screen->get_dims(1),
            (int)d_phase->d_screen->get_dims(2),
            (int)ps->d_tscreen->d_screen->get_dims(1),
            xoff[std::make_pair("atmos", idx)],
            yoff[std::make_pair("atmos", idx)], this->block_size);
      } else {
        if (this->lgs) {
          delta = 1.0f - ps->altitude / this->d_lgs->hg;
        }
        target_raytrace(this->d_phase->d_screen->get_data(),
                          ps->d_tscreen->d_screen->get_data(),
                          (int)d_phase->d_screen->get_dims(1),
                          (int)d_phase->d_screen->get_dims(2),
                          (int)ps->d_tscreen->d_screen->get_dims(1),
                          xoff[std::make_pair("atmos", idx)] - ps->accumx,
                          yoff[std::make_pair("atmos", idx)] - ps->accumy, this->G,
                          this->thetaML, this->dx, this->dy, this->block_size, delta);
      }
    }
    p++;
  }

  return EXIT_SUCCESS;
}

int SutraSource::raytrace(SutraDms *ydms, bool rst, bool do_phase_var,
                           bool async) {
  current_context->set_active_device(device, 1);
  if (rst) this->d_phase->d_screen->reset();
  map<type_screen, float>::iterator p;
  p = xoff.begin();
  float delta = 1.0f;
  while (p != xoff.end()) {
    string types = p->first.first;
    if ((types.find("pzt") == 0) || (types.find("tt") == 0) ||
        (types.find("kl") == 0)) {
      int inddm = p->first.second;
      if (inddm < 0) throw "error in SutraSource::raytrace, dm not find";
      SutraDm *ps = ydms->d_dms[inddm];
      if ((p == xoff.end()) && async) {
        target_raytrace_async(this->phase_telemetry,
                              this->d_phase->d_screen->get_data(),
                              ps->d_shape->d_screen->get_data(),
                              (int)d_phase->d_screen->get_dims(1),
                              (int)d_phase->d_screen->get_dims(2),
                              (int)ps->d_shape->d_screen->get_dims(1),
                              xoff[make_pair(types, inddm)],
                              yoff[make_pair(types, inddm)], this->block_size);
      } else {
        if (this->lgs) {
          delta = 1.0f - ps->altitude / this->d_lgs->hg;
        }
        target_raytrace(this->d_phase->d_screen->get_data(),
                          ps->d_shape->d_screen->get_data(),
                          (int)d_phase->d_screen->get_dims(1),
                          (int)d_phase->d_screen->get_dims(2),
                          (int)ps->d_shape->d_screen->get_dims(1),
                          xoff[std::make_pair(types, inddm)],
                          yoff[std::make_pair(types, inddm)], this->G * ps->G,
                          this->thetaML + ps->thetaML, this->dx + ps->dx, this->dy + ps->dy, this->block_size, delta);
        }
      }
    p++;
  }

  if (type != "wfs") {
    // select phase pixels in the valid portion of pupil
    fillindx(this->d_phasepts->get_data(), this->d_phase->d_screen->get_data(),
             this->d_wherephase->get_data(), this->scale,
             this->d_wherephase->get_nb_elements(),
             current_context->get_device(device));

    float phase_avg = 0;
    // compute avg phase in the pupil
    phase_avg = this->d_phasepts->sum();
    phase_avg /= this->d_wherephase->get_nb_elements();

    // substract avg from phase in the pupil
    fillindx(this->d_phasepts->get_data(), this->d_phase->d_screen->get_data(),
             this->d_wherephase->get_data(), this->scale, -phase_avg,
             this->d_wherephase->get_nb_elements(),
             current_context->get_device(device));

    // compute instantaneous phase variance and average
    this->phase_var = this->d_phasepts->dot(this->d_phasepts, 1, 1);
    this->phase_var /= this->d_wherephase->get_nb_elements();
    if (do_phase_var) {
      if (this->phase_var_count >= 0) this->phase_var_avg += this->phase_var;

      this->phase_var_count += 1;
    }
  }

  return EXIT_SUCCESS;
}

int SutraSource::raytrace(bool rst) {
  this->current_context->set_active_device(this->device, 1);

  if (rst) {
    this->d_phase->d_screen->reset();
  }

  if (this->d_ncpa_phase != nullptr) {
    this->d_phase->d_screen->axpy(1.0f, this->d_ncpa_phase, 1, 1);
  }
  return EXIT_SUCCESS;
}

int SutraSource::raytrace(SutraTelescope *tel, bool rst) {
  this->current_context->set_active_device(this->device, 1);

  if (rst) {
    this->d_phase->d_screen->reset();
  }

  if (tel->d_phase_ab_M1 != nullptr && this->type != "wfs") {
    this->d_phase->d_screen->axpy(1.0f, tel->d_phase_ab_M1, 1, 1);
  }
  if (tel->d_phase_ab_M1_m != nullptr && this->type == "wfs") {
    this->d_phase->d_screen->axpy(1.0f, tel->d_phase_ab_M1_m, 1, 1);
  }
  return EXIT_SUCCESS;
}

int SutraSource::raytrace(SutraTelescope *tel, SutraAtmos *atmos,
                           SutraDms *ydms, bool do_phase_var, bool async) {
  this->current_context->set_active_device(this->device, 1);
  this->raytrace(atmos, async);
  this->raytrace(tel, false);
  this->raytrace(ydms, false, do_phase_var, async);
  this->raytrace(false);
  return EXIT_SUCCESS;
}
int SutraSource::comp_image(int puponly, bool comp_le) {
  current_context->set_active_device(device, 1);
  if (this->d_amplipup == 0) return -1;

  // set complex amplitude in the pupil plane to zero
  carma_safe_call(
      cudaMemset(this->d_amplipup->get_data(), 0,
                 sizeof(cuFloatComplex) * this->d_amplipup->get_nb_elements()));

  // fill complex amplitude in the pupil with phase @ lambda
  fill_amplipup(
      this->d_amplipup->get_data(), this->d_phase->d_screen->get_data(),
      this->d_pupil->get_data(), this->scale, puponly,
      this->d_phase->d_screen->get_dims(1), this->d_phase->d_screen->get_dims(2),
      this->d_amplipup->get_dims(1), current_context->get_device(device));

  // do fft : complex amplitude in the focal plane
  CarmaFFT(this->d_amplipup->get_data(), this->d_amplipup->get_data(), 1,
            *this->d_amplipup->get_plan());

  // take square norm to retreive short exposure image
  abs2(this->d_image_se->get_data(), this->d_amplipup->get_data(),
       this->d_image_se->get_dims(1) * this->d_image_se->get_dims(2),
       current_context->get_device(device));

  // scale image because of fft
  this->d_image_se->scale(1.0f / this->d_wherephase->get_dims(1), 1);
  this->d_image_se->scale(1.0f / this->d_wherephase->get_dims(1), 1);

  // if long exposure image is not null
  if (this->d_image_le != 0L && comp_le) {
    // add new short exposure
    this->d_image_le->axpy(1.0f, this->d_image_se, 1, 1);
    this->strehl_counter += 1;
  }
  return EXIT_SUCCESS;
}

float sinc(double x) {
  if (x == 0) return 1;
  x *= CARMA_PI;
  return sin(x) / x;
}

/**
 * @brief fit the strehl with a sinc
 *
 * Utilise la “croix” de 3 pixels centraux qui encadrent le max
 * pour fitter des paraboles qui determinent la position du maximum,
 * puis calcule l’interpolation exacte en ce point via la formule
 * des sinus cardinaux qui s’applique a un signal bien echantillonne.
 *
 * @param d_img full image of size img_size*img_size
 * @param ind_max position of the maximum in d_img
 * @param img_size size of the d_img leading dimension
 * @return float Strehl fitted
 */
float SutraSource::fitmax2x1dSinc(float *d_img, int ind_max, int img_size) {
  const int small_size = 2 * this->d_smallimg_size + 1;
  extract(this->d_smallimg->get_data(), d_img, img_size, ind_max, small_size,
          true);
  std::vector<float> smallimg(small_size * small_size);
  const int center = ((small_size * small_size) - 1) / 2;
  const int right_ind = center + 1;
  const int left_ind = center - 1;
  const int down_ind = center - small_size;
  const int up_ind = center + small_size;

  this->d_smallimg->device2host(smallimg.data());

  const float A =
      0.5 * (smallimg[down_ind] + smallimg[up_ind]) - smallimg[center];
  const float B =
      0.5 * (smallimg[left_ind] + smallimg[right_ind]) - smallimg[center];
  if (A * B == 0) {
    DEBUG_TRACE("ERROR: can not estimate the SR");
    return 0;
  }
  const float x0 = -0.5 * (smallimg[up_ind] - smallimg[down_ind]) / 2. / A;
  const float y0 = -0.5 * (smallimg[right_ind] - smallimg[left_ind]) / 2. / B;

  // interpolation exacte par shannon(theoreme des sinus cardinaux)
  float valmax = 0.0;
  int ind = 0;
  for (int i = -this->d_smallimg_size; i < this->d_smallimg_size + 1; ++i) {
    const float tmpi = sinc(x0 - i);
    for (int j = -this->d_smallimg_size; j < this->d_smallimg_size + 1; ++j) {
      const float tmpj = sinc(y0 - j);
      valmax += tmpi * tmpj * smallimg[ind++];
    }
  }
  return valmax;
}

int SutraSource::comp_strehl(bool do_fit) {
  current_context->set_active_device(device, 1);

  const int max_se = this->d_image_se->aimax(1);
  const int max_le = this->d_image_le->aimax(1);
  if (do_fit) {
    const int img_size = d_image_se->get_dims(1);
    this->strehl_se =
        fitmax2x1dSinc(this->d_image_se->get_data(), max_se, img_size);
    this->strehl_le =
        fitmax2x1dSinc(this->d_image_le->get_data(), max_le, img_size);
  } else {
    cudaMemcpy(&(this->strehl_se), &(this->d_image_se->get_data()[max_se]),
               sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(&(this->strehl_le), &(this->d_image_le->get_data()[max_le]),
               sizeof(float), cudaMemcpyDeviceToHost);
  }

  if (this->strehl_counter > 0)
    this->strehl_le /= this->strehl_counter;
  else
    this->strehl_le = this->strehl_se;

  return EXIT_SUCCESS;
}
