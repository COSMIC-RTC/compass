#include <sutra_rtc.h>
#include <sutra_wfs_geom.h>
#include <sutra_wfs_sh.h>
#include <sutra_wfs_pyr_pyrhr.h>

sutra_rtc::sutra_rtc(carma_context *context) {
  this->current_context = context;
  this->device = context->get_activeDevice();
}

sutra_rtc::~sutra_rtc() {
  current_context->set_activeDevice(device,1);

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

  //delete this->current_context;
}
int sutra_rtc::add_centroider(int nwfs, long nvalid, float offset, float scale, long device, char *typec) {
  sutra_sensors *sensors = nullptr;
  add_centroider(sensors, nwfs, nvalid, offset, scale, device, typec);

  return EXIT_SUCCESS;
}

int sutra_rtc::add_centroider(sutra_sensors *sensors, int nwfs, long nvalid, float offset, float scale,
                              long device, char *typec) {
  current_context->set_activeDevice(device,1);
  if (strcmp(typec, "bpcog") == 0)
    d_centro.push_back(
      new sutra_centroider_bpcog(current_context, sensors, nwfs, nvalid, offset, scale,
                                 device, 10));
  else if (strcmp(typec, "cog") == 0)
    d_centro.push_back(
      new sutra_centroider_cog(current_context, sensors, nwfs, nvalid, offset, scale,
                               device));
  else if (strcmp(typec, "corr") == 0)
    d_centro.push_back(
      new sutra_centroider_corr(current_context, sensors, nwfs, nvalid, offset, scale,
                                device));
  else if (strcmp(typec, "pyr") == 0)
    d_centro.push_back(
      new sutra_centroider_pyr(current_context, sensors, nwfs, nvalid, offset, scale,
                               device));
  else if (strcmp(typec, "roof") == 0)
    d_centro.push_back(
      new sutra_centroider_roof(current_context, sensors, nwfs, nvalid, offset, scale,
                                device));
  else if (strcmp(typec, "tcog") == 0)
    d_centro.push_back(
      new sutra_centroider_tcog(current_context, sensors, nwfs, nvalid, offset, scale,
                                device));
  else if (strcmp(typec, "wcog") == 0)
    d_centro.push_back(
      new sutra_centroider_wcog(current_context, sensors, nwfs, nvalid, offset, scale,
                                device));
  else {
    std::cerr << "centroider unknown" << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

int sutra_rtc::add_controller_geo(int nactu, int Nphi, float delay,
                                  long device, sutra_dms *dms, char **type_dmseen, float *alt, int ndm, bool wfs_direction) {
  current_context->set_activeDevice(device,1);
  this->d_control.push_back(
    new sutra_controller_geo(current_context, nactu, Nphi, delay, dms, type_dmseen, alt, ndm, wfs_direction));
  return EXIT_SUCCESS;
}

int sutra_rtc::add_controller(int nactu, float delay, long device, const char *typec) {
  current_context->set_activeDevice(device,1);
  int ncentroids = 0;
  for (size_t idx = 0; idx < (this->d_centro).size(); idx++)
    ncentroids += this->d_centro[idx]->nvalid;

  string type_ctr(typec);
  sutra_dms *dms = nullptr;
  char **type_dmseen = nullptr;
  float *alt = nullptr;

  if (type_ctr.compare("generic") == 0) {
    d_control.push_back(
      new sutra_controller_generic(current_context, ncentroids, nactu, delay, dms, type_dmseen, alt, 0));
  }

  return EXIT_SUCCESS;
}
int sutra_rtc::add_controller(int nactu, float delay, long device,
                              const char *typec, sutra_dms *dms, char **type_dmseen, float *alt, int ndm) {
  current_context->set_activeDevice(device,1);
  int ncentroids = 0;
  for (size_t idx = 0; idx < (this->d_centro).size(); idx++)
    ncentroids += this->d_centro[idx]->nvalid;

  current_context->set_activeDevice(device,1);
  string type_ctr(typec);
  if (type_ctr.compare("ls") == 0) {
    d_control.push_back(
      new sutra_controller_ls(current_context, ncentroids, nactu, delay, dms, type_dmseen, alt, ndm));
  } else if (type_ctr.compare("cured") == 0) {
    d_control.push_back(
      new sutra_controller_cured(current_context, ncentroids, nactu, delay, dms, type_dmseen, alt, ndm));
  } else if (type_ctr.compare("mv") == 0) {
    d_control.push_back(
      new sutra_controller_mv(current_context, ncentroids, nactu, delay, dms, type_dmseen, alt, ndm));
  } else if (type_ctr.compare("generic") == 0) {
    d_control.push_back(
      new sutra_controller_generic(current_context, ncentroids, nactu, delay, dms, type_dmseen, alt, ndm));
  } else if ((type_ctr.compare("kalman_GPU") == 0)
             || (type_ctr.compare("kalman_CPU") == 0)) {
    d_control.push_back(
      new sutra_controller_kalman(current_context, ncentroids, nactu, dms, type_dmseen, alt, ndm));
  } else {
    DEBUG_TRACE("Controller '%s' unknown\n", typec);
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

int sutra_rtc::rm_controller() {

  current_context->set_activeDevice(device,1);
  //for (size_t idx = 0; idx < (this->d_control).size(); idx++) {
  while ((this->d_control).size() > 0) {
    delete this->d_control.back();
    d_control.pop_back();
  }

  return EXIT_SUCCESS;
}

int sutra_rtc::do_imat(int ncntrl, sutra_dms *ydm) {
  current_context->set_activeDevice(device,1);
  carma_obj<float> *d_imat = NULL;
  if (this->d_control[ncntrl]->get_type().compare("ls") == 0) {
    SCAST(sutra_controller_ls *, control, this->d_control[ncntrl]);
    d_imat = control->d_imat;
  } else if (this->d_control[ncntrl]->get_type().compare("mv") == 0) {
    SCAST(sutra_controller_mv *, control, this->d_control[ncntrl]);
    d_imat = control->d_imat;
  } else {
    DEBUG_TRACE("Controller needs to be ls or mv\n");
    return EXIT_SUCCESS;
  }

  vector<sutra_dm *>::iterator p;
  p = this->d_control[ncntrl]->d_dmseen.begin();
  int inds1 = 0;
  int cc = 0;
  int cc2 = 0;

  std::cout << "Doing imat..." << std::endl;
  while (p != this->d_control[ncntrl]->d_dmseen.end()) {
    sutra_dm *dm = *p;
    std::string desc = "DM"+carma_utils::to_string(cc2);
    auto progress = carma_utils::ProgressBar(dm->ninflu, desc);
    for (int j = 0; j < dm->ninflu; ++j) {
      // Push
      dm->comp_oneactu(j, dm->push4imat);

      for (size_t idx_cntr = 0; idx_cntr < (this->d_centro).size();
           ++idx_cntr) {
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
      //DEBUG_TRACE("actu %d", j);
      do_centroids(ncntrl, true);

      int device = this->d_control[ncntrl]->d_centroids->getDevice();
      convert_centro(*this->d_control[ncntrl]->d_centroids,
                     *this->d_control[ncntrl]->d_centroids, 0, 0.5f / dm->push4imat,
                     this->d_control[ncntrl]->d_centroids->getNbElem(),
                     current_context->get_device(device));
      this->d_control[ncntrl]->d_centroids->copyInto(d_imat->getDataAt(inds1),
          this->d_control[ncntrl]->nslope());
      dm->reset_shape();
      // Pull
      dm->comp_oneactu(j, -1.0f * dm->push4imat);
      for (size_t idx_cntr = 0; idx_cntr < (this->d_centro).size();
           idx_cntr++) {
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
      device = this->d_control[ncntrl]->d_centroids->getDevice();
      do_centroids(ncntrl, true);
      convert_centro(*this->d_control[ncntrl]->d_centroids,
                     *this->d_control[ncntrl]->d_centroids, 0, 0.5f / dm->push4imat,
                     this->d_control[ncntrl]->d_centroids->getNbElem(),
                     current_context->get_device(device));

      float alphai = -1.0f;
      cublasSaxpy(current_context->get_cublasHandle(),
                  this->d_control[ncntrl]->d_centroids->getNbElem(), &alphai,
                  this->d_control[ncntrl]->d_centroids->getData(), 1, d_imat->getDataAt(inds1),
                  1);

      dm->reset_shape();
      inds1 += this->d_control[ncntrl]->nslope();
      ++cc;
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

int sutra_rtc::do_imat_geom(int ncntrl, sutra_dms *ydm,
                            int type) {
  current_context->set_activeDevice(device,1);
  vector<sutra_dm *>::iterator p;
  p = this->d_control[ncntrl]->d_dmseen.begin();
  int inds1, inds2;
  inds1 = 0;
  while (p != this->d_control[ncntrl]->d_dmseen.end()) {
    sutra_dm *dm = *p;
    for (int j = 0; j < dm->ninflu; j++) {
      dm->comp_oneactu(j, dm->push4imat); //
      inds2 = 0;
      for (size_t idx_cntr = 0; idx_cntr < (this->d_control).size();
           idx_cntr++) {

        carma_obj<float> *d_imat;
        if (this->d_control[ncntrl]->get_type().compare("ls") == 0) {
          sutra_controller_ls * control = dynamic_cast<sutra_controller_ls *>(this->d_control[ncntrl]);
          d_imat = control->d_imat;
        } else if (this->d_control[ncntrl]->get_type().compare("mv") == 0) {
          sutra_controller_mv * control = dynamic_cast<sutra_controller_mv *>(this->d_control[ncntrl]);
          d_imat = control->d_imat;
        } else {
          DEBUG_TRACE("controller should be ls or mv");
          return EXIT_FAILURE;
        }

        sutra_wfs *wfs = this->d_centro[idx_cntr]->wfs;

        wfs->sensor_trace(ydm, 1);
        //sensors->d_wfs[nwfs]->comp_image();
        if(wfs->type != "sh") {
          sutra_wfs_sh *_wfs = dynamic_cast<sutra_wfs_sh *>(wfs);
          _wfs->slopes_geom(type, d_imat[inds1 + inds2]);
        } else if(wfs->type != "geo") {
          sutra_wfs_geom *_wfs = dynamic_cast<sutra_wfs_geom *>(wfs);
          _wfs->slopes_geom(type, d_imat[inds1 + inds2]);
        } else if(wfs->type != "pyrhr") {
          sutra_wfs_pyr_pyrhr *_wfs = dynamic_cast<sutra_wfs_pyr_pyrhr *>(wfs);
          _wfs->slopes_geom(type, d_imat[inds1 + inds2]);
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

int sutra_rtc::remove_ref(int ncntrl) {
  carma_axpy<float>(this->current_context->get_cublasHandle(),
                    this->d_control[ncntrl]->d_centroids->getNbElem(),
                    -1.0f,
                    this->d_control[ncntrl]->d_centroids_ref->getData(), 1,
                    this->d_control[ncntrl]->d_centroids->getData(), 1);
  return EXIT_SUCCESS;
}

int sutra_rtc::do_centroids() {
  current_context->set_activeDevice(device,1);
  for (size_t idx_cntr = 0; idx_cntr < (this->d_centro).size(); idx_cntr++) {
    this->d_centro[idx_cntr]->get_cog();
  }

  return EXIT_SUCCESS;
}

int sutra_rtc::do_centroids(int ncntrl) {
  return do_centroids(ncntrl, true);
}

int sutra_rtc::do_centroids(int ncntrl, bool noise) {
  current_context->set_activeDevice(device,1);
  int indssp = 0;

  for (size_t idx_cntr = 0; idx_cntr < (this->d_centro).size(); idx_cntr++) {

    this->d_centro[idx_cntr]->get_cog(this->d_control[ncntrl]->d_subsum->getDataAt(indssp),
                                      this->d_control[ncntrl]->d_centroids->getDataAt(2*indssp),noise);

    indssp += this->d_centro[idx_cntr]->wfs->nvalid_tot;
  }
  remove_ref(ncntrl);

  return EXIT_SUCCESS;
}

int sutra_rtc::do_centroids(int ncntrl, float *bincube, int npix, int ntot) {
  carma_streams *streams = nullptr;
  int indssp = 0;

  for (size_t idx_cntr = 0; idx_cntr < (this->d_centro).size(); idx_cntr++) {
    this->d_centro[idx_cntr]->get_cog(streams, bincube, this->d_control[ncntrl]->d_subsum->getDataAt(indssp),
                                      this->d_control[ncntrl]->d_centroids->getDataAt(2*indssp),this->d_centro[idx_cntr]->nvalid, npix, ntot);

    indssp += this->d_centro[idx_cntr]->nvalid;
  }
  remove_ref(ncntrl);

  return EXIT_SUCCESS;
}
int sutra_rtc::do_centroids_geom(int ncntrl) {
  current_context->set_activeDevice(device,1);
  int inds2 = 0;

  for (size_t idx_cntr = 0; idx_cntr < (this->d_centro).size(); idx_cntr++) {

    sutra_wfs *wfs = this->d_centro[idx_cntr]->wfs;
    if(wfs->type == "sh") {
      sutra_wfs_sh *_wfs = dynamic_cast<sutra_wfs_sh *>(wfs);
      _wfs->slopes_geom(0,this->d_control[ncntrl]->d_centroids->getDataAt(inds2));

    } else if(wfs->type == "geo") {
      sutra_wfs_geom *_wfs = dynamic_cast<sutra_wfs_geom *>(wfs);
      _wfs->slopes_geom(0,this->d_control[ncntrl]->d_centroids->getDataAt(inds2));

    } else if(wfs->type == "pyrhr") {
      sutra_wfs_pyr_pyrhr *_wfs = dynamic_cast<sutra_wfs_pyr_pyrhr *>(wfs);
      _wfs->slopes_geom(0,this->d_control[ncntrl]->d_centroids->getDataAt(inds2));
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

int sutra_rtc::set_centroids_ref(int ncntrl, float *centroids_ref) {
  current_context->set_activeDevice(device,1);
  this->d_control[ncntrl]->set_centroids_ref(centroids_ref);
  return EXIT_SUCCESS;
}

int sutra_rtc::get_centroids_ref(int ncntrl, float *centroids_ref) {
  current_context->set_activeDevice(device,1);
  this->d_control[ncntrl]->get_centroids_ref(centroids_ref);
  return EXIT_SUCCESS;
}

int sutra_rtc::do_control(int ncntrl) {
  current_context->set_activeDevice(device,1);

  if(this->d_control[ncntrl]->open_loop) {
    return EXIT_SUCCESS;
  } else
    this->d_control[ncntrl]->comp_com();

  return EXIT_SUCCESS;
}

int sutra_rtc::do_clipping(int ncntrl, float min, float max) {
  current_context->set_activeDevice(device,1);
  this->d_control[ncntrl]->clip_voltage(min, max);

  return EXIT_SUCCESS;
}

int sutra_rtc::comp_voltage(int ncntrl) {
  this->d_control[ncntrl]->comp_voltage();
  this->d_control[ncntrl]->command_delay();

  return EXIT_SUCCESS;
}
int sutra_rtc::apply_control(int ncntrl, sutra_dms *ydm) {
  current_context->set_activeDevice(device,1);

  comp_voltage(ncntrl);

  vector<sutra_dm *>::iterator p;
  p = this->d_control[ncntrl]->d_dmseen.begin();
  int idx = 0;
  if ((this->d_control[ncntrl]->get_type().compare("ls") == 0)
      || (this->d_control[ncntrl]->get_type().compare("mv") == 0)
      || (this->d_control[ncntrl]->get_type().compare("geo") == 0)) {
    // "streamed" controllers case

    while (p != this->d_control[ncntrl]->d_dmseen.end()) {
      sutra_dm *dm = *p;

      int nstreams = this->d_control[ncntrl]->streams->get_nbStreams();
      if (nstreams > dm->ninflu) {
        for (int i = 0; i < nstreams; i++) {
          int istart = i * dm->ninflu / nstreams;
          carmaSafeCall(
            cudaMemcpyAsync(dm->d_comm->getDataAt(istart),
                            this->d_control[ncntrl]->d_voltage->getDataAt(idx + istart),
                            sizeof(float) * dm->ninflu / nstreams,
                            cudaMemcpyDeviceToDevice,
                            (*this->d_control[ncntrl]->streams)[i]));
          dm->comp_shape();
        }
      } else {
        dm->comp_shape(this->d_control[ncntrl]->d_voltage->getDataAt(idx));
      }
      idx += dm->ninflu;
      p++;
    }
  } else { // "non-streamed" controllers
    while (p != this->d_control[ncntrl]->d_dmseen.end()) {
      sutra_dm *dm = *p;
      dm->comp_shape(this->d_control[ncntrl]->d_voltage->getDataAt(idx));
      idx += dm->ninflu;
      p++;
    }
  }
  return EXIT_SUCCESS;
}
