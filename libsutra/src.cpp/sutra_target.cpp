#include <sutra_ao_utils.h>
#include <sutra_target.h>
#include <sutra_phase.h>

int fft_goodsize(long size) {
  int mradix = 2;
  float tmpf = 0.;
  long tmpl = 0;

  tmpf = logf(size) / logf(3);
  tmpl = (long) tmpf;
  tmpf -= tmpl;
  mradix = (
      tmpf > (logf(size) / logf(mradix) - (long) (logf(size) / logf(mradix))) ?
          3 : mradix);

  tmpf = logf(size) / logf(5);
  tmpl = (long) tmpf;
  tmpf -= tmpl;
  mradix = (
      tmpf > (logf(size) / logf(mradix) - (long) (logf(size) / logf(mradix))) ?
          5 : mradix);

  tmpf = logf(size) / logf(7);
  tmpl = (long) tmpf;
  tmpf -= tmpl;
  mradix = (
      tmpf > (logf(size) / logf(mradix) - (long) (logf(size) / logf(mradix))) ?
          7 : mradix);

  return mradix;
}

sutra_source::sutra_source(carma_context *context, float xpos, float ypos,
    float lambda, float mag, long size, string type, float *pupil, int device) {

  this->init_source(context, xpos, ypos, lambda, mag, size, type, device);
  this->d_pupil = new carma_obj<float>(this->current_context,
      this->d_phase->d_screen->getDims());
  this->d_pupil->host2device(pupil);
  long sumpup = (long) this->d_pupil->asum(1);
  long *dims_data1 = new long[2];
  dims_data1[0] = 1;
  dims_data1[1] = pow(2, (long) (logf(sumpup) / logf(2)) + 1);
  this->d_phasepts = new carma_obj<float>(this->current_context, dims_data1);
  dims_data1[1] = sumpup;
  this->d_wherephase = new carma_obj<int>(this->current_context, dims_data1);
  int *wherephase = new int[sumpup];
  int cpt = 0;
  for (int cc = 0; cc < this->d_pupil->getNbElem(); cc++) {
    if (pupil[cc] > 0) {
      wherephase[cpt] = cc;
      cpt += 1;
    }
  }

  this->d_wherephase->host2device(wherephase);
  delete[] wherephase;
  delete[] dims_data1;
}

sutra_source::sutra_source(carma_context *context, float xpos, float ypos,
    float lambda, float mag, long size, string type, int device) {
  this->init_source(context, xpos, ypos, lambda, mag, size, type, device);
}

inline int sutra_source::init_source(carma_context *context, float xpos,
    float ypos, float lambda, float mag, long size, string type, int device) {
  this->current_context = context;
  this->strehl_counter = 0;

  this->tposx = xpos;
  this->tposy = ypos;
  this->lambda = lambda;
  this->mag = mag;
  this->npos = size;

  this->d_phase = new sutra_phase(context, size);
    
  // ADDING INSTRUMENTAL PHASE
  // this->d_phase_instru = new sutra_phase(context, size);
  //

  this->type = type;
  this->device = device;
  this->scale = float(2 * 3.14159265 / lambda); // phase is expected in microns

  //struct cudaDeviceProp deviceProperties;
  //cudaGetDeviceProperties(&deviceProperties, device);
  //this->blockSize = (int)sqrt(deviceProperties.maxThreadsPerBlock);
  this->block_size = 8;
  this->phase_var_avg = 0;
  this->phase_var_count = 0;

  long *dims_data2 = new long[3];
  dims_data2[0] = 2;
  dims_data2[1] = size;
  dims_data2[2] = size;
  int nstreams = (size + this->block_size - size % this->block_size)
      / this->block_size;

  this->phase_telemetry = new carma_host_obj<float>(dims_data2, MA_WRICOMB,
      nstreams);

  this->d_image = 0L;
  this->d_pupil = 0L;
  this->d_amplipup = 0L;
  this->d_leimage = 0L;
  this->d_phasepts = 0L;
  this->d_wherephase = 0L;

  if (type != "wfs") {
    int mradix = 2; //fft_goodsize(size);

    int fft_size = pow(mradix, (long) (logf(2 * size) / logf(mradix)) + 1);
    dims_data2[1] = fft_size;
    dims_data2[2] = fft_size;

    this->d_image = new carma_obj<float>(context, dims_data2);
    this->d_amplipup = new carma_obj<cuFloatComplex>(context, dims_data2);

    cufftHandle *plan = this->d_amplipup->getPlan(); ///< FFT plan
    cufftSafeCall(
        cufftPlan2d(plan, this->d_amplipup->getDims(1),
            this->d_amplipup->getDims(2), CUFFT_C2C));
  }

  this->lgs = false;

  delete[] dims_data2;

  return EXIT_SUCCESS;

}

sutra_source::~sutra_source() {
  //delete this->current_context;

  delete this->d_phase;
    
  // REMOVING PHASE INSTRU
  // delete this->d_phase_instru;
  //
    
  delete this->phase_telemetry;

  if (d_image != 0L)
    delete this->d_image;
  if (d_pupil != 0L)
    delete this->d_pupil;
  if (d_leimage != 0L)
    delete this->d_leimage;
  if (d_amplipup != 0L)
    delete this->d_amplipup;
  if (d_phasepts != 0L)
    delete this->d_phasepts;
  if (d_wherephase != 0L)
    delete this->d_wherephase;
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

int sutra_source::init_strehlmeter() {
  this->strehl_counter = 0;
  this->comp_image(1);

  cudaMemcpy(&(this->ref_strehl),
      &(this->d_image->getData()[this->d_image->imax(1) - 1]), sizeof(float),
      cudaMemcpyDeviceToHost);

  if (this->d_leimage == 0L)
    this->d_leimage = new carma_obj<float>(this->current_context,
        this->d_image->getDims());
  else
    cutilSafeCall(
        cudaMemset(this->d_leimage->getData(), 0,
            sizeof(float) * this->d_leimage->getNbElem()));

  return EXIT_SUCCESS;
}

int sutra_source::add_layer(string type, float alt, float mxoff, float myoff) {
  xoff[make_pair(type, alt)] = mxoff;
  yoff[make_pair(type, alt)] = myoff;

  return EXIT_SUCCESS;
}

int sutra_source::remove_layer(string type, float alt) {
  xoff.erase(make_pair(type, alt));
  yoff.erase(make_pair(type, alt));

  return EXIT_SUCCESS;
}

int sutra_source::raytrace_shm(sutra_atmos *yatmos) {
  cutilSafeCall(
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
            (int) d_phase->d_screen->getDims(1),
            (int) d_phase->d_screen->getDims(2),
            (int) ps->second->d_tscreen->d_screen->getDims(1),
            (int) ps->second->d_tscreen->d_screen->getDims(2),
            xoff[make_pair("atmos", alt)], yoff[make_pair("atmos", alt)],
            ps->second->d_tscreen->d_screen->getNbElem(),
            ps->second->channelDesc, current_context->get_device(device));
      }
    }
    p++;
  }

  return EXIT_SUCCESS;
}

int sutra_source::raytrace(sutra_atmos *yatmos, bool async) {
//  cutilSafeCall(cudaDeviceSynchronize());
  cutilSafeCall(
      cudaMemset(this->d_phase->d_screen->getData(), 0,
          sizeof(float) * this->d_phase->d_screen->getNbElem()));

  float delta;
  map<type_screen, float>::iterator p;
  p = xoff.begin();

  while (p != xoff.end()) {
    string types = p->first.first;

    if (types.find("atmos") == 0) {
      float alt = p->first.second;
      sutra_tscreen * ps;
      ps = yatmos->d_screens[alt];
      p++;
      if ((p == xoff.end()) && async) {

        target_raytrace_async(this->phase_telemetry,
            this->d_phase->d_screen->getData(),
            ps->d_tscreen->d_screen->getData(),
            (int) d_phase->d_screen->getDims(1),
            (int) d_phase->d_screen->getDims(2),
            (int) ps->d_tscreen->d_screen->getDims(1),
            xoff[make_pair("atmos", alt)], yoff[make_pair("atmos", alt)],
            this->block_size);
      } else {
    	  if (this->lgs){
    		  delta = 1.0f - alt/this->d_lgs->hg;
    		  target_lgs_raytrace(this->d_phase->d_screen->getData(),
    		              ps->d_tscreen->d_screen->getData(),
    		              (int) d_phase->d_screen->getDims(1),
    		              (int) d_phase->d_screen->getDims(2),
    		              (int) ps->d_tscreen->d_screen->getDims(1),
    		              xoff[make_pair(types, alt)], yoff[make_pair(types, alt)],delta,
    		              this->block_size);
    	  }
    	  else

    		  target_raytrace(this->d_phase->d_screen->getData(),
    				  ps->d_tscreen->d_screen->getData(),
    				  (int) d_phase->d_screen->getDims(1),
    				  (int) d_phase->d_screen->getDims(2),
    				  (int) ps->d_tscreen->d_screen->getDims(1),
    				  xoff[make_pair("atmos", alt)], yoff[make_pair("atmos", alt)],
    				  this->block_size);
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

int sutra_source::raytrace(sutra_dms *ydms, int rst, bool async) {
  if (rst == 1)
    cutilSafeCall(
        cudaMemset(this->d_phase->d_screen->getData(), 0,
            sizeof(float) * this->d_phase->d_screen->getNbElem()));
  map<type_screen, float>::iterator p;
  p = xoff.begin();
  while (p != xoff.end()) {
    string types = p->first.first;
    if ((types.find("pzt") == 0) || (types.find("tt") == 0)
        || (types.find("kl") == 0)) {
      float alt = p->first.second;
      sutra_dm *ps = ydms->d_dms[make_pair(types, alt)];
      p++;
      if ((p == xoff.end()) && async) {
        target_raytrace_async(this->phase_telemetry,
            this->d_phase->d_screen->getData(),
            ps->d_shape->d_screen->getData(),
            (int) d_phase->d_screen->getDims(1),
            (int) d_phase->d_screen->getDims(2),
            (int) ps->d_shape->d_screen->getDims(1),
            xoff[make_pair(types, alt)], yoff[make_pair(types, alt)],
            this->block_size);
      } else{
    		  target_raytrace(this->d_phase->d_screen->getData(),
    				  ps->d_shape->d_screen->getData(),
    				  (int) d_phase->d_screen->getDims(1),
    				  (int) d_phase->d_screen->getDims(2),
    				  (int) ps->d_shape->d_screen->getDims(1),
    				  xoff[make_pair(types, alt)], yoff[make_pair(types, alt)],
    				  this->block_size);
      }
    } else
      p++;
  }

  if (type != "wfs") {
    // select phase pixels in the valid portion of pupil
    fillindx(this->d_phasepts->getData(), this->d_phase->d_screen->getData(),
        this->d_wherephase->getData(), this->scale,
        this->d_wherephase->getNbElem(), current_context->get_device(device));

    float phase_avg = 0;
    // compute avg phase in the pupil
    phase_avg = this->d_phasepts->sum();
    phase_avg /= this->d_wherephase->getNbElem();

    // substract avg from phase in the pupil
    fillindx(this->d_phasepts->getData(), this->d_phase->d_screen->getData(),
        this->d_wherephase->getData(), this->scale, -phase_avg,
        this->d_wherephase->getNbElem(), current_context->get_device(device));

    // compute instantaneous phase variance and average
    this->phase_var = this->d_phasepts->dot(this->d_phasepts, 1, 1);
    this->phase_var /= this->d_wherephase->getNbElem();

    this->phase_var_avg += this->phase_var;

    this->phase_var_count += 1;
  }

  return EXIT_SUCCESS;
}

int sutra_source::raytrace(sutra_dms *ydms, int rst) {
  raytrace(ydms, rst, false);

  return EXIT_SUCCESS;
}

int sutra_source::comp_image(int puponly) {
  if (this->d_amplipup == 0)
    return -1;

  // set complex amplitude in the pupil plane to zero
  cutilSafeCall(
      cudaMemset(this->d_amplipup->getData(), 0,
          sizeof(cuFloatComplex) * this->d_amplipup->getNbElem()));

  /*
   fillpupil(this->d_amplipup->getData(), mask, this->d_phase->d_screen->getDims(1),
   this->d_phase->d_screen->getDims(2), this->d_amplipup->getDims(1),current_context->get_device(device));

   */

  // fill complex amplitude in the pupil with phase @ lambda
  fill_amplipup(this->d_amplipup->getData(), this->d_phase->d_screen->getData(),
      this->d_pupil->getData(), this->scale, puponly,
      this->d_phase->d_screen->getDims(1), this->d_phase->d_screen->getDims(2),
      this->d_amplipup->getDims(1), current_context->get_device(device));
    
//  // fill complex amplitude in the pupil with phase @ lambda
//  fill_amplipup(this->d_amplipup->getData(), this->d_phase->d_screen->getData(),
//              this->d_phase_instru->d_screen_instru->getData(),
//              this->d_pupil->getData(), this->scale, puponly,
//              this->d_phase->d_screen->getDims(1), this->d_phase->d_screen->getDims(2),
//              this->d_amplipup->getDims(1), current_context->get_device(device));

  // do fft : complex amplitude in the focal plane
  carma_fft(this->d_amplipup->getData(), this->d_amplipup->getData(), 1,
      *this->d_amplipup->getPlan());

  // take square norm to retreive short exposure image
  abs2(this->d_image->getData(), this->d_amplipup->getData(),
      this->d_image->getDims(1) * this->d_image->getDims(2), current_context->get_device(device));

  // scale image because of fft
  //this->d_image->scale(1.0f / this->d_image->getNbElem(), 1);

  // if long exposure image is not null
  if (this->d_leimage != 0L) {
    // add new short exposure
    this->d_leimage->axpy(1.0f, this->d_image, 1, 1);
    this->strehl_counter += 1;
  }

  return EXIT_SUCCESS;
}

int sutra_source::comp_strehl() {
  //this->strehl_se = expf(-this->phase_var);
  //this->strehl_le = expf(-this->phase_var_avg/this->phase_var_count);
  
  cudaMemcpy(&(this->strehl_se),
        &(this->d_image->getData()[this->d_image->imax(1) - 1]), sizeof(float),
        cudaMemcpyDeviceToHost);
  cudaMemcpy(&(this->strehl_le),
        &(this->d_leimage->getData()[this->d_leimage->imax(1) - 1]),
        sizeof(float), cudaMemcpyDeviceToHost);
  this->strehl_se /= (this->d_wherephase->getDims(1) * this->d_wherephase->getDims(1));
  this->strehl_le /=  (this->strehl_counter * ( this->d_wherephase->getDims(1)* this->d_wherephase->getDims(1)));
  
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

sutra_target::sutra_target(carma_context *context, int ntargets, float *xpos,
    float *ypos, float *lambda, float *mag, long *sizes, float *pupil,
    int device) {
  this->ntargets = ntargets;

  for (int i = 0; i < ntargets; i++) {
    d_targets.push_back(
        new sutra_source(context, xpos[i], ypos[i], lambda[i], mag[i], sizes[i],
            "target", pupil, device));
  }
}

sutra_target::~sutra_target() {
//  for (size_t idx = 0; idx < (this->d_targets).size(); idx++) {
    while((this->d_targets).size()>0) {
    delete this->d_targets.back();
    d_targets.pop_back();
  }
}

