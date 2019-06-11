#include <sutra_source.h>

sutra_source::sutra_source(carma_context *context, float xpos, float ypos,
                           float lambda, float mag, float zerop, long size,
                           string type, carma_obj<float> *pupil, int Npts,
                           int device)
    : current_context(context),
      device(device),
      d_pupil(pupil),
      d_ncpa_phase(nullptr) {
  current_context->set_activeDevice(device, 1);

  this->init_source(context, xpos, ypos, lambda, mag, zerop, size, type,
                    device);
  // float h_pupil[this->d_pupil->getNbElem()];
  float *h_pupil = new float[this->d_pupil->getNbElem()];
  this->d_pupil->device2host(h_pupil);

  long *dims_data1 = new long[2];
  dims_data1[0] = 1;
  dims_data1[1] = pow(2, (long)(logf(Npts) / logf(2)) + 1);
  // this->d_phasepts = new carma_obj<float>(this->current_context, dims_data1);
  dims_data1[1] = Npts;
  this->d_phasepts = new carma_obj<float>(this->current_context, dims_data1);
  this->d_wherephase = new carma_obj<int>(this->current_context, dims_data1);
  int *wherephase = new int[Npts];
  int cpt = 0;
  for (int cc = 0; cc < this->d_pupil->getNbElem(); cc++) {
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

sutra_source::sutra_source(carma_context *context, float xpos, float ypos,
                           float lambda, float mag, float zerop, long size,
                           string type, int device)
    : current_context(context), device(device), d_ncpa_phase(nullptr) {
  this->device = device;
  current_context = context;
  current_context->set_activeDevice(device, 1);
  this->init_source(context, xpos, ypos, lambda, mag, zerop, size, type,
                    device);
}

inline int sutra_source::init_source(carma_context *context, float xpos,
                                     float ypos, float lambda, float mag,
                                     float zerop, long size, string type,
                                     int device) {
  current_context->set_activeDevice(device, 1);
  this->current_context = context;
  this->strehl_counter = 0;

  this->posx = xpos;
  this->posy = ypos;
  this->lambda = lambda;
  this->mag = mag;
  this->zp = zerop;
  this->npts = size;

  this->d_phase = new sutra_phase(context, size);

  // ADDING INSTRUMENTAL PHASE
  // this->d_phase_instru = new sutra_phase(context, size);
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
      new carma_host_obj<float>(dims_data2, MA_WRICOMB, nstreams);

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

    this->d_image_se = new carma_obj<float>(context, dims_data2);
    this->d_amplipup = new carma_obj<cuFloatComplex>(context, dims_data2);

    cufftHandle *plan = this->d_amplipup->getPlan();  ///< FFT plan
    carmafftSafeCall(cufftPlan2d(plan, this->d_amplipup->getDims(1),
                                 this->d_amplipup->getDims(2), CUFFT_C2C));
  }

  this->lgs = false;
  this->G = 1.0f;
  this->thetaML = 0.0f;
  this->dx = 0.0f;
  this->dy = 0.0f;

  long *dims_data1 = new long[2];
  dims_data1[0] = 1;

  /// temporary array for accurate strehl computation
  dims_data1[1] = 9;
  this->d_smallimg = new carma_obj<float>(context, dims_data1);

  /// temporary array for accurate strehl computation
  dims_data2[1] = 9;
  dims_data2[2] = 6;
  std::vector<float> fitmat = {
      // v
      18,  -36, 18,  18,  -36, 18,  18,  -36, 18,  18,  18,  18, -36, -36,
      -36, 18,  18,  18,  27,  0,   -27, 0,   0,   0,   -27, 0,  27,  -18,
      0,   18,  -18, 0,   18,  -18, 0,   18,  -18, -18, -18, 0,  0,   0,
      18,  18,  18,  -12, 24,  -12, 24,  60,  24,  -12, 24,  -12};
  // std::vector<float> fitmat = { // v.T
  //     18,  18, 27,  -18, -18, -12, -36, 18,  0,  0,  -18, 24,  18, 18,
  //     -27, 18, -18, -12, 18,  -36, 0,   -18, 0,  24, -36, -36, 0,  0,
  //     0,   60, 18,  -36, 0,   18,  0,   24,  18, 18, -27, -18, 18, -12,
  //     -36, 18, 0,   0,   18,  24,  18,  18,  27, 18, 18,  -12};
  this->d_fitmat = new carma_obj<float>(context, dims_data2, fitmat.data());

  /// temporary array for accurate strehl computation
  dims_data1[1] = 6;
  this->d_fitpoly = new carma_obj<float>(context, dims_data1);

  delete[] dims_data2;
  delete[] dims_data1;

  return EXIT_SUCCESS;
}

sutra_source::~sutra_source() {
  // delete this->current_context;
  this->current_context->set_activeDevice(this->device, 1);
  delete this->d_smallimg;
  delete this->d_fitmat;
  delete this->d_fitpoly;
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

int sutra_source::init_strehlmeter() {
  current_context->set_activeDevice(device, 1);
  this->strehl_counter = 0;
  this->comp_image(1, false);

  cudaMemcpy(&(this->ref_strehl),
             &(this->d_image_se->getData()[this->d_image_se->aimax(1)]),
             sizeof(float), cudaMemcpyDeviceToHost);

  if (this->d_image_le == 0L)
    this->d_image_le = new carma_obj<float>(this->current_context,
                                            this->d_image_se->getDims());
  else {  // Reset strehl case
    carmaSafeCall(cudaMemset(this->d_image_le->getData(), 0,
                             sizeof(float) * this->d_image_le->getNbElem()));
    this->strehl_counter = 0;
    this->phase_var_avg = 0.f;
    this->phase_var_count = 0;
  }

  return EXIT_SUCCESS;
}

int sutra_source::reset_strehlmeter() {
  carmaSafeCall(cudaMemset(this->d_image_le->getData(), 0,
                           sizeof(float) * this->d_image_le->getNbElem()));
  this->strehl_counter = 0;

  return EXIT_SUCCESS;
}

int sutra_source::reset_phase() {
  this->d_phase->d_screen->reset();
  return EXIT_SUCCESS;
}
int sutra_source::add_layer(string type, int idx, float mxoff, float myoff) {
  current_context->set_activeDevice(device, 1);
  xoff[make_pair(type, idx)] = mxoff;
  yoff[make_pair(type, idx)] = myoff;

  return EXIT_SUCCESS;
}

int sutra_source::remove_layer(string type, int idx) {
  current_context->set_activeDevice(device, 1);
  xoff.erase(make_pair(type, idx));
  yoff.erase(make_pair(type, idx));

  return EXIT_SUCCESS;
}

int sutra_source::raytrace(sutra_atmos *yatmos, bool async) {
  //  carmaSafeCall(cudaDeviceSynchronize());
  current_context->set_activeDevice(device, 1);
  carmaSafeCall(
      cudaMemset(this->d_phase->d_screen->getData(), 0,
                 sizeof(float) * this->d_phase->d_screen->getNbElem()));

  float delta;
  map<type_screen, float>::iterator p;
  p = xoff.begin();

  while (p != xoff.end()) {
    string types = p->first.first;

    if (types.find("atmos") == 0) {
      int idx = p->first.second;
      sutra_tscreen *ps;
      ps = yatmos->d_screens[idx];
      p++;
      if ((p == xoff.end()) && async) {
        target_raytrace_async(
            this->phase_telemetry, this->d_phase->d_screen->getData(),
            ps->d_tscreen->d_screen->getData(),
            (int)d_phase->d_screen->getDims(1),
            (int)d_phase->d_screen->getDims(2),
            (int)ps->d_tscreen->d_screen->getDims(1),
            xoff[std::make_pair("atmos", idx)],
            yoff[std::make_pair("atmos", idx)], this->block_size);
      } else {
        if (this->lgs) {
          delta = 1.0f - ps->altitude / this->d_lgs->hg;
          if (delta > 0)
            target_lgs_raytrace(this->d_phase->d_screen->getData(),
                                ps->d_tscreen->d_screen->getData(),
                                (int)d_phase->d_screen->getDims(1),
                                (int)d_phase->d_screen->getDims(2),
                                (int)ps->d_tscreen->d_screen->getDims(1),
                                xoff[make_pair(types, idx)],
                                yoff[make_pair(types, idx)], delta,
                                this->block_size);
        } else

          target_raytrace(this->d_phase->d_screen->getData(),
                          ps->d_tscreen->d_screen->getData(),
                          (int)d_phase->d_screen->getDims(1),
                          (int)d_phase->d_screen->getDims(2),
                          (int)ps->d_tscreen->d_screen->getDims(1),
                          xoff[std::make_pair("atmos", idx)],
                          yoff[std::make_pair("atmos", idx)], this->G,
                          this->thetaML, this->dx, this->dy, this->block_size);
      }
    } else
      p++;
  }

  return EXIT_SUCCESS;
}

int sutra_source::raytrace(sutra_dms *ydms, bool rst, bool do_phase_var,
                           bool async) {
  current_context->set_activeDevice(device, 1);
  if (rst) this->d_phase->d_screen->reset();
  map<type_screen, float>::iterator p;
  p = xoff.begin();
  while (p != xoff.end()) {
    string types = p->first.first;
    if ((types.find("pzt") == 0) || (types.find("tt") == 0) ||
        (types.find("kl") == 0)) {
      int inddm = p->first.second;
      if (inddm < 0) throw "error in sutra_source::raytrace, dm not find";
      sutra_dm *ps = ydms->d_dms[inddm];
      p++;
      if ((p == xoff.end()) && async) {
        target_raytrace_async(this->phase_telemetry,
                              this->d_phase->d_screen->getData(),
                              ps->d_shape->d_screen->getData(),
                              (int)d_phase->d_screen->getDims(1),
                              (int)d_phase->d_screen->getDims(2),
                              (int)ps->d_shape->d_screen->getDims(1),
                              xoff[make_pair(types, inddm)],
                              yoff[make_pair(types, inddm)], this->block_size);
      } else {
        target_raytrace(this->d_phase->d_screen->getData(),
                        ps->d_shape->d_screen->getData(),
                        (int)d_phase->d_screen->getDims(1),
                        (int)d_phase->d_screen->getDims(2),
                        (int)ps->d_shape->d_screen->getDims(1),
                        xoff[std::make_pair(types, inddm)],
                        yoff[std::make_pair(types, inddm)], this->G,
                        this->thetaML, this->dx, this->dy, this->block_size);
      }
    } else
      p++;
  }

  if (type != "wfs") {
    // select phase pixels in the valid portion of pupil
    fillindx(this->d_phasepts->getData(), this->d_phase->d_screen->getData(),
             this->d_wherephase->getData(), this->scale,
             this->d_wherephase->getNbElem(),
             current_context->get_device(device));

    float phase_avg = 0;
    // compute avg phase in the pupil
    phase_avg = this->d_phasepts->sum();
    phase_avg /= this->d_wherephase->getNbElem();

    // substract avg from phase in the pupil
    fillindx(this->d_phasepts->getData(), this->d_phase->d_screen->getData(),
             this->d_wherephase->getData(), this->scale, -phase_avg,
             this->d_wherephase->getNbElem(),
             current_context->get_device(device));

    // compute instantaneous phase variance and average
    this->phase_var = this->d_phasepts->dot(this->d_phasepts, 1, 1);
    this->phase_var /= this->d_wherephase->getNbElem();
    if (do_phase_var) {
      if (this->phase_var_count >= 0) this->phase_var_avg += this->phase_var;

      this->phase_var_count += 1;
    }
  }

  return EXIT_SUCCESS;
}

int sutra_source::raytrace(bool rst) {
  this->current_context->set_activeDevice(this->device, 1);

  if (rst) {
    this->d_phase->d_screen->reset();
  }

  if (this->d_ncpa_phase != nullptr) {
    this->d_phase->d_screen->axpy(1.0f, this->d_ncpa_phase, 1, 1);
  }
  return EXIT_SUCCESS;
}

int sutra_source::raytrace(sutra_telescope *tel, bool rst) {
  this->current_context->set_activeDevice(this->device, 1);

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

int sutra_source::raytrace(sutra_telescope *tel, sutra_atmos *atmos,
                           sutra_dms *ydms, bool do_phase_var, bool async) {
  this->current_context->set_activeDevice(this->device, 1);
  this->raytrace(atmos, async);
  this->raytrace(tel, false);
  this->raytrace(ydms, false, do_phase_var, async);
  this->raytrace(false);
  return EXIT_SUCCESS;
}
int sutra_source::comp_image(int puponly, bool comp_le) {
  current_context->set_activeDevice(device, 1);
  if (this->d_amplipup == 0) return -1;

  // set complex amplitude in the pupil plane to zero
  carmaSafeCall(
      cudaMemset(this->d_amplipup->getData(), 0,
                 sizeof(cuFloatComplex) * this->d_amplipup->getNbElem()));

  // fill complex amplitude in the pupil with phase @ lambda
  fill_amplipup(
      this->d_amplipup->getData(), this->d_phase->d_screen->getData(),
      this->d_pupil->getData(), this->scale, puponly,
      this->d_phase->d_screen->getDims(1), this->d_phase->d_screen->getDims(2),
      this->d_amplipup->getDims(1), current_context->get_device(device));

  // do fft : complex amplitude in the focal plane
  carma_fft(this->d_amplipup->getData(), this->d_amplipup->getData(), 1,
            *this->d_amplipup->getPlan());

  // take square norm to retreive short exposure image
  abs2(this->d_image_se->getData(), this->d_amplipup->getData(),
       this->d_image_se->getDims(1) * this->d_image_se->getDims(2),
       current_context->get_device(device));

  // scale image because of fft
  this->d_image_se->scale(1.0f / this->d_wherephase->getDims(1), 1);
  this->d_image_se->scale(1.0f / this->d_wherephase->getDims(1), 1);

  // if long exposure image is not null
  if (this->d_image_le != 0L && comp_le) {
    // add new short exposure
    this->d_image_le->axpy(1.0f, this->d_image_se, 1, 1);
    this->strehl_counter += 1;
  }
  return EXIT_SUCCESS;
}

float sutra_source::fitmax(float *d_img, int ind_max, int img_size) {
  // int max_y = ind_max / img_size;
  // int max_x = ind_max % img_size;
  // DEBUG_TRACE("ind_max %d, max_x %d, max_y %d", ind_max, max_x, max_y);

  extract(this->d_smallimg->getData(), d_img, img_size, ind_max, 3, true);

  // DEBUG_TRACE("smallimg_krnl");
  // std::vector<float> smallimg(3 * 3);
  // this->d_smallimg->device2host(smallimg.data());
  // for (int i = 0; i < 9; ++i) {
  //   std::cerr << smallimg[i] << ", ";
  // }
  // std::cerr << "\n";

  // int index = 0;
  // for (int j = -1; j < 2; ++j) {
  //   int y = max_y + j;
  //   if (y < 0) y += img_size;
  //   if (y > img_size - 1) y -= img_size;
  //   for (int i = -1; i < 2; ++i) {
  //     int x = max_x + i;
  //     if (x < 0) x += img_size;
  //     if (x > img_size - 1) x -= img_size;
  //     cudaMemcpy(this->d_smallimg->getData() + index,
  //                &(d_img[x + y * img_size]), sizeof(float),
  //                cudaMemcpyDeviceToDevice);
  //     ++index;
  //   }
  // }

  // DEBUG_TRACE("smallimg_cpy");
  // std::vector<float> smallimg2(3 * 3);
  // this->d_smallimg->device2host(smallimg2.data());
  // for (int i = 0; i < 9; ++i) {
  //   std::cerr << smallimg2[i] << ", ";
  // }
  // std::cerr << "\n";

  // DEBUG_TRACE("fitmat");
  // std::vector<float> fitmat(6*9);
  // this->d_fitmat->device2host(fitmat.data());
  // for(int i=0; i<6*9; ++i){
  //   std::cerr << fitmat[i] << " ";
  // }
  // std::cerr << "\n";

  this->d_fitpoly->gemv('t', 1. / 108., this->d_fitmat, 9, this->d_smallimg, 1,
                        0., 1);
  float fitpoly[6];
  this->d_fitpoly->device2host(fitpoly);
  const float A = fitpoly[0];
  const float B = fitpoly[1];
  const float C = fitpoly[2];
  const float denom = C * C - 4 * A * B;
  float x0 = 0.0, y0 = 0.0;
  if (denom != 0) {
    y0 = (2. * B * fitpoly[3] - fitpoly[4] * C) / denom;
    x0 = (2. * A * fitpoly[4] - fitpoly[3] * C) / denom;
  }

  // DEBUG_TRACE("fitpoly");
  // for (int i = 0; i < 6; ++i) {
  //   std::cerr << fitpoly[i] << ", ";
  // }
  // std::cerr << "\n";

  float D = fitpoly[5] - (x0 * y0 * C + A * x0 * x0 + B * y0 * y0);
  // DEBUG_TRACE("%f %f %f", x0, y0, D);
  return D;
}

int sutra_source::comp_strehl(bool do_fit) {
  current_context->set_activeDevice(device, 1);

  const int max_se = this->d_image_se->aimax(1);
  const int max_le = this->d_image_le->aimax(1);
  if (do_fit) {
    const int img_size = d_image_se->getDims(1);
    this->strehl_se = fitmax(this->d_image_se->getData(), max_se, img_size);
    this->strehl_le = fitmax(this->d_image_le->getData(), max_le, img_size);
  } else {
    cudaMemcpy(&(this->strehl_se), &(this->d_image_se->getData()[max_se]),
               sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(&(this->strehl_le), &(this->d_image_le->getData()[max_le]),
               sizeof(float), cudaMemcpyDeviceToHost);
  }

  // DEBUG_TRACE(
  //     "%f %f %f %f", this->strehl_se, this->strehl_le,
  //     fitmax(this->d_image_se->getData(), max_se, d_image_se->getDims(1)),
  //     fitmax(this->d_image_le->getData(), max_le, d_image_le->getDims(1)));

  // this->strehl_se /= (this->d_wherephase->getDims(1) *
  // this->d_wherephase->getDims(1));
  if (this->strehl_counter > 0)
    this->strehl_le /= this->strehl_counter;
  else
    this->strehl_le = this->strehl_se;

  /*
  this->strehl_se /= this->ref_strehl;
  this->strehl_le /= (this->ref_strehl * this->strehl_counter);*/

  return EXIT_SUCCESS;
}
