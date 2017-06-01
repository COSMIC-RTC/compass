#include <sutra_ao_utils.h>
#include <sutra_phase.h>
#include <sutra_target.h>
#include <sutra_telescope.h>

int fft_goodsize(long size) {
  int mradix = 2;
  float tmpf = 0.;
  long tmpl = 0;

  tmpf = logf(size) / logf(3);
  tmpl = (long)tmpf;
  tmpf -= tmpl;
  mradix =
      (tmpf > (logf(size) / logf(mradix) - (long)(logf(size) / logf(mradix)))
           ? 3
           : mradix);

  tmpf = logf(size) / logf(5);
  tmpl = (long)tmpf;
  tmpf -= tmpl;
  mradix =
      (tmpf > (logf(size) / logf(mradix) - (long)(logf(size) / logf(mradix)))
           ? 5
           : mradix);

  tmpf = logf(size) / logf(7);
  tmpl = (long)tmpf;
  tmpf -= tmpl;
  mradix =
      (tmpf > (logf(size) / logf(mradix) - (long)(logf(size) / logf(mradix)))
           ? 7
           : mradix);

  return mradix;
}

sutra_source::sutra_source(carma_context *context, float xpos, float ypos,
                           float lambda, float mag, float zerop, long size,
                           string type, carma_obj<float> *pupil, int Npts,
                           int device)
    : current_context(context), device(device), d_pupil(pupil),
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

  this->tposx = xpos;
  this->tposy = ypos;
  this->lambda = lambda;
  this->mag = mag;
  this->zp = zerop;
  this->npos = size;

  this->d_phase = new sutra_phase(context, size);

  // ADDING INSTRUMENTAL PHASE
  // this->d_phase_instru = new sutra_phase(context, size);
  //

  this->type = type;
  this->device = device;
  this->scale = float(2 * 3.14159265 / lambda); // phase is expected in microns

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

  this->d_image = 0L;
  this->d_amplipup = 0L;
  this->d_leimage = 0L;
  this->d_phasepts = 0L;
  this->d_wherephase = 0L;

  if (type != "wfs") {
    int mradix = 2; // fft_goodsize(size);

    int fft_size = pow(mradix, (long)(logf(2 * size) / logf(mradix)) + 1);
    dims_data2[1] = fft_size;
    dims_data2[2] = fft_size;

    this->d_image = new carma_obj<float>(context, dims_data2);
    this->d_amplipup = new carma_obj<cuFloatComplex>(context, dims_data2);

    cufftHandle *plan = this->d_amplipup->getPlan(); ///< FFT plan
    carmafftSafeCall(cufftPlan2d(plan, this->d_amplipup->getDims(1),
                                 this->d_amplipup->getDims(2), CUFFT_C2C));
  }

  this->lgs = false;
  this->G = 1.0f;
  this->thetaML = 0.0f;
  this->dx = 0.0f;
  this->dy = 0.0f;

  delete[] dims_data2;

  return EXIT_SUCCESS;
}

sutra_source::~sutra_source() {
  // delete this->current_context;
  this->current_context->set_activeDevice(this->device, 1);

  delete this->d_phase;

  // REMOVING PHASE INSTRU
  // delete this->d_phase_instru;
  //

  delete this->phase_telemetry;

  if (this->d_image != 0L)
    delete this->d_image;
  if (this->d_leimage != 0L)
    delete this->d_leimage;
  if (this->d_amplipup != 0L)
    delete this->d_amplipup;
  if (this->d_phasepts != 0L)
    delete this->d_phasepts;
  if (this->d_wherephase != 0L)
    delete this->d_wherephase;
  if (this->d_ncpa_phase != nullptr)
    delete this->d_ncpa_phase;
  /*
   for( std::map< type_screen,float>::iterator it = this->xoff.begin();
   this->xoff.end()!= it ; it++) {
   delete it->second;
   }
   for( std::map< type_screen,float>::iterator it = this->yoff.begin();
   this->yoff.end()!= it ; it++) {
   delete it->second;
   }
   */

  this->xoff.clear();
  this->yoff.clear();
}

int sutra_source::get_ncpa_phase(float *h_src, size_t size) {
  size_t phase_size = this->d_phase->d_screen->getNbElem();
  if (size != phase_size) {
    DEBUG_TRACE("ERROR: argument given is incorrect, wrong size (%ld given, "
                "%ld expected)",
                size, phase_size);
    return EXIT_FAILURE;
  }
  if (this->d_ncpa_phase == nullptr) {
    DEBUG_TRACE("ERROR: no NCPA has been loaded");
    return EXIT_FAILURE;
  }
  this->d_ncpa_phase->device2host(h_src);
  return EXIT_SUCCESS;
}

int sutra_source::set_ncpa_phase(float *h_dest, size_t size) {
  size_t phase_size = this->d_phase->d_screen->getNbElem();
  if (size != phase_size) {
    DEBUG_TRACE("ERROR: argument given is incorrect, wrong size (%ld given, "
                "%ld expected)",
                size, phase_size);
    return EXIT_FAILURE;
  }
  if (this->d_ncpa_phase == nullptr) {
    this->d_ncpa_phase = new carma_obj<float>(this->d_phase->d_screen);
  }
  this->d_ncpa_phase->host2device(h_dest);
  return EXIT_SUCCESS;
}

int sutra_source::init_strehlmeter() {
  current_context->set_activeDevice(device, 1);
  this->strehl_counter = 0;
  this->comp_image(1, false);

  cudaMemcpy(&(this->ref_strehl),
             &(this->d_image->getData()[this->d_image->imax(1) - 1]),
             sizeof(float), cudaMemcpyDeviceToHost);

  if (this->d_leimage == 0L)
    this->d_leimage =
        new carma_obj<float>(this->current_context, this->d_image->getDims());
  else { // Reset strehl case
    carmaSafeCall(cudaMemset(this->d_leimage->getData(), 0,
                             sizeof(float) * this->d_leimage->getNbElem()));
    this->strehl_counter = 0;
    this->phase_var_avg = 0.f;
    this->phase_var_count = 0;
  }

  return EXIT_SUCCESS;
}

int sutra_source::reset_strehlmeter() {

  carmaSafeCall(cudaMemset(this->d_leimage->getData(), 0,
                           sizeof(float) * this->d_leimage->getNbElem()));
  this->strehl_counter = 0;

  return EXIT_SUCCESS;
}

int sutra_source::add_layer(string type, float alt, float mxoff, float myoff) {
  current_context->set_activeDevice(device, 1);
  xoff[make_pair(type, alt)] = mxoff;
  yoff[make_pair(type, alt)] = myoff;

  return EXIT_SUCCESS;
}

int sutra_source::remove_layer(string type, float alt) {
  current_context->set_activeDevice(device, 1);
  xoff.erase(make_pair(type, alt));
  yoff.erase(make_pair(type, alt));

  return EXIT_SUCCESS;
}

int sutra_source::raytrace_shm(sutra_atmos *yatmos) {
  current_context->set_activeDevice(device, 1);
  carmaSafeCall(
      cudaMemset(this->d_phase->d_screen->getData(), 0,
                 sizeof(float) * this->d_phase->d_screen->getNbElem()));

  map<type_screen, float>::iterator p;
  p = xoff.begin();
  while (p != xoff.end()) {
    string types = p->first.first;
    if (types.find("atmos") == 0) {
      float alt = p->first.second;
      map<float, sutra_tscreen *>::iterator ps;
      ps = yatmos->d_screens.find(alt);
      if (ps != yatmos->d_screens.end()) {
        target_texraytrace(this->d_phase->d_screen->getData(),
                           ps->second->d_tscreen->d_screen->getData(),
                           (int)d_phase->d_screen->getDims(1),
                           (int)d_phase->d_screen->getDims(2),
                           (int)ps->second->d_tscreen->d_screen->getDims(1),
                           (int)ps->second->d_tscreen->d_screen->getDims(2),
                           xoff[std::make_pair("atmos", alt)],
                           yoff[std::make_pair("atmos", alt)],
                           ps->second->d_tscreen->d_screen->getNbElem(),
                           ps->second->channelDesc,
                           current_context->get_device(device));
      }
    }
    p++;
  }

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
      float alt = p->first.second;
      sutra_tscreen *ps;
      ps = yatmos->d_screens[alt];
      p++;
      if ((p == xoff.end()) && async) {

        target_raytrace_async(
            this->phase_telemetry, this->d_phase->d_screen->getData(),
            ps->d_tscreen->d_screen->getData(),
            (int)d_phase->d_screen->getDims(1),
            (int)d_phase->d_screen->getDims(2),
            (int)ps->d_tscreen->d_screen->getDims(1),
            xoff[std::make_pair("atmos", alt)],
            yoff[std::make_pair("atmos", alt)], this->block_size);
      } else {
        if (this->lgs) {
          delta = 1.0f - alt / this->d_lgs->hg;
          if (delta > 0)
            target_lgs_raytrace(this->d_phase->d_screen->getData(),
                                ps->d_tscreen->d_screen->getData(),
                                (int)d_phase->d_screen->getDims(1),
                                (int)d_phase->d_screen->getDims(2),
                                (int)ps->d_tscreen->d_screen->getDims(1),
                                xoff[make_pair(types, alt)],
                                yoff[make_pair(types, alt)], delta,
                                this->block_size);
        } else

          target_raytrace(this->d_phase->d_screen->getData(),
                          ps->d_tscreen->d_screen->getData(),
                          (int)d_phase->d_screen->getDims(1),
                          (int)d_phase->d_screen->getDims(2),
                          (int)ps->d_tscreen->d_screen->getDims(1),
                          xoff[std::make_pair("atmos", alt)],
                          yoff[std::make_pair("atmos", alt)], this->G,
                          this->thetaML, this->dx, this->dy, this->block_size);
      }
    } else
      p++;
  }

  return EXIT_SUCCESS;
}

int sutra_source::raytrace(sutra_atmos *yatmos) {
  raytrace(yatmos, false);

  return EXIT_SUCCESS;
}

int sutra_source::raytrace(sutra_dms *ydms, int rst, bool async,
                           int do_phase_var) {
  current_context->set_activeDevice(device, 1);
  if (rst == 1)
    this->d_phase->d_screen->reset();
  map<type_screen, float>::iterator p;
  p = xoff.begin();
  while (p != xoff.end()) {
    string types = p->first.first;
    if ((types.find("pzt") == 0) || (types.find("tt") == 0) ||
        (types.find("kl") == 0)) {
      float alt = p->first.second;
      int inddm = ydms->get_inddm(types, alt);
      if (inddm < 0)
        throw "error in sutra_source::raytrace, dm not find";
      sutra_dm *ps = ydms->d_dms[inddm];
      p++;
      if ((p == xoff.end()) && async) {
        target_raytrace_async(
            this->phase_telemetry, this->d_phase->d_screen->getData(),
            ps->d_shape->d_screen->getData(),
            (int)d_phase->d_screen->getDims(1),
            (int)d_phase->d_screen->getDims(2),
            (int)ps->d_shape->d_screen->getDims(1), xoff[make_pair(types, alt)],
            yoff[make_pair(types, alt)], this->block_size);
      } else {
        target_raytrace(this->d_phase->d_screen->getData(),
                        ps->d_shape->d_screen->getData(),
                        (int)d_phase->d_screen->getDims(1),
                        (int)d_phase->d_screen->getDims(2),
                        (int)ps->d_shape->d_screen->getDims(1),
                        xoff[std::make_pair(types, alt)],
                        yoff[std::make_pair(types, alt)], this->G,
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
      if (this->phase_var_count >= 0)
        this->phase_var_avg += this->phase_var;

      this->phase_var_count += 1;
    }
  }

  return EXIT_SUCCESS;
}

int sutra_source::raytrace(sutra_dms *ydms, int rst, int phase_var) {
  raytrace(ydms, rst, false, phase_var);
  return EXIT_SUCCESS;
}

int sutra_source::raytrace(int rst) {
  if (this->d_ncpa_phase == nullptr) {
    // DEBUG_TRACE("ERROR: no NCPA has been loaded");
    return EXIT_FAILURE;
  }

  this->current_context->set_activeDevice(this->device, 1);

  // if (rst == 1)
  //   this->d_phase->d_screen->reset();
  //
  // target_raytrace(
  //     this->d_phase->d_screen->getData(), this->d_ncpa_phase->getData(),
  //     this->d_phase->d_screen->getDims(1), this->d_phase->d_screen->getDims(2),
  //     this->d_ncpa_phase->getDims(1), 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f,
  //     this->block_size);

  float alpha = (rst == 1)?0.f:1.f;
  this->d_phase->d_screen->axpy(alpha, this->d_ncpa_phase, 1, 1);

  return EXIT_SUCCESS;
}

int sutra_source::comp_image(int puponly, bool comp_le) {
  current_context->set_activeDevice(device, 1);
  if (this->d_amplipup == 0)
    return -1;

  // set complex amplitude in the pupil plane to zero
  carmaSafeCall(
      cudaMemset(this->d_amplipup->getData(), 0,
                 sizeof(cuFloatComplex) * this->d_amplipup->getNbElem()));

  /*
   fillpupil(this->d_amplipup->getData(), mask,
   this->d_phase->d_screen->getDims(1),
   this->d_phase->d_screen->getDims(2),
   this->d_amplipup->getDims(1),current_context->get_device(device));

   */

  // fill complex amplitude in the pupil with phase @ lambda
  fill_amplipup(
      this->d_amplipup->getData(), this->d_phase->d_screen->getData(),
      this->d_pupil->getData(), this->scale, puponly,
      this->d_phase->d_screen->getDims(1), this->d_phase->d_screen->getDims(2),
      this->d_amplipup->getDims(1), current_context->get_device(device));

  //  // fill complex amplitude in the pupil with phase @ lambda
  //  fill_amplipup(this->d_amplipup->getData(),
  //  this->d_phase->d_screen->getData(),
  //              this->d_phase_instru->d_screen_instru->getData(),
  //              this->d_pupil->getData(), this->scale, puponly,
  //              this->d_phase->d_screen->getDims(1),
  //              this->d_phase->d_screen->getDims(2),
  //              this->d_amplipup->getDims(1),
  //              current_context->get_device(device));

  // do fft : complex amplitude in the focal plane
  carma_fft(this->d_amplipup->getData(), this->d_amplipup->getData(), 1,
            *this->d_amplipup->getPlan());

  // take square norm to retreive short exposure image
  abs2(this->d_image->getData(), this->d_amplipup->getData(),
       this->d_image->getDims(1) * this->d_image->getDims(2),
       current_context->get_device(device));

  // scale image because of fft
  this->d_image->scale(1.0f / this->d_wherephase->getDims(1), 1);
  this->d_image->scale(1.0f / this->d_wherephase->getDims(1), 1);

  // if long exposure image is not null
  if (this->d_leimage != 0L && comp_le) {
    // add new short exposure
    this->d_leimage->axpy(1.0f, this->d_image, 1, 1);
    this->strehl_counter += 1;
  }
  return EXIT_SUCCESS;
}

int sutra_source::comp_strehl() {
  // this->strehl_se = expf(-this->phase_var);
  // this->strehl_le = expf(-this->phase_var_avg/this->phase_var_count);

  current_context->set_activeDevice(device, 1);
  cudaMemcpy(&(this->strehl_se),
             &(this->d_image->getData()[this->d_image->imax(1) - 1]),
             sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(&(this->strehl_le),
             &(this->d_leimage->getData()[this->d_leimage->imax(1) - 1]),
             sizeof(float), cudaMemcpyDeviceToHost);
  // this->strehl_se /= (this->d_wherephase->getDims(1) *
  // this->d_wherephase->getDims(1));
  if (this->strehl_counter > 0)
    this->strehl_le /= this->strehl_counter;
  else
    this->strehl_le = this->strehl_se;

  /*
  cudaMemcpy(&(this->strehl_se),
      &(this->d_image->getData()[this->d_image->imax(1) - 1]), sizeof(float),
      cudaMemcpyDeviceToHost);
  cudaMemcpy(&(this->strehl_le),
      &(this->d_leimage->getData()[this->d_leimage->imax(1) - 1]),
      sizeof(float), cudaMemcpyDeviceToHost);

  this->strehl_se /= this->ref_strehl;
  this->strehl_le /= (this->ref_strehl * this->strehl_counter);*/

  return EXIT_SUCCESS;
}

sutra_target::sutra_target(carma_context *context, sutra_telescope *yTelescope,
                           int ntargets, float *xpos, float *ypos,
                           float *lambda, float *mag, float zerop, long *sizes,
                           int Npts, int device) {

  this->ntargets = ntargets;

  for (int i = 0; i < ntargets; i++) {
    d_targets.push_back(new sutra_source(context, xpos[i], ypos[i], lambda[i],
                                         mag[i], zerop, sizes[i], "target",
                                         yTelescope->d_pupil, Npts, device));
  }
}

sutra_target::~sutra_target() {
  //  for (size_t idx = 0; idx < (this->d_targets).size(); idx++) {
  while ((this->d_targets).size() > 0) {
    delete this->d_targets.back();
    d_targets.pop_back();
  }
}
