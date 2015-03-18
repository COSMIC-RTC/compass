#include <sutra_rtc.h>
#include <sutra_wfs_geom.h>
#include <sutra_wfs_sh.h>

sutra_rtc::sutra_rtc(carma_context *context) {
  this->current_context = context;
  this->device = context->get_activeDevice();
  this->delay = 0;
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
    cerr << "centroider unknown\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

int sutra_rtc::add_controller_geo(int nactu, int Nphi, long delay,
    long device, sutra_dms *dms, char **type_dmseen, float *alt, int ndm) {
  current_context->set_activeDevice(device,1);
  this->d_control.push_back(
      new sutra_controller_geo(current_context, nactu, Nphi, delay, dms, type_dmseen, alt, ndm));
  return EXIT_SUCCESS;
}

int sutra_rtc::add_controller(int nactu, long delay, long device,
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

  map<type_screen, sutra_dm *>::iterator p;
  p = this->d_control[ncntrl]->d_dmseen.begin();
  int inds1;
  inds1 = 0;
  while (p != this->d_control[ncntrl]->d_dmseen.end()) {
    for (int j = 0; j < p->second->ninflu; j++) {
      // Push
      p->second->comp_oneactu(j, p->second->push4imat);

      for (size_t idx_cntr = 0; idx_cntr < (this->d_centro).size();
          idx_cntr++) {
        sutra_wfs *wfs = this->d_centro[idx_cntr]->wfs;
        float tmp_noise = wfs->noise;
        wfs->noise = -1;
        wfs->kernconv = true;
        wfs->sensor_trace(ydm, 1);
        wfs->comp_image();
        wfs->noise = tmp_noise;
        wfs->kernconv = false;
      }
      //DEBUG_TRACE("actu %d", j);
      do_centroids(ncntrl, true);

      int device = this->d_control[ncntrl]->d_centroids->getDevice();
      convert_centro(*this->d_control[ncntrl]->d_centroids,
          *this->d_control[ncntrl]->d_centroids, 0, 0.5f / p->second->push4imat,
          this->d_control[ncntrl]->d_centroids->getNbElem(),
          current_context->get_device(device));
      this->d_control[ncntrl]->d_centroids->copyInto((*d_imat)[inds1],
          this->d_control[ncntrl]->nslope());
      p->second->reset_shape();
      // Pull
      p->second->comp_oneactu(j, -1.0f * p->second->push4imat);
      for (size_t idx_cntr = 0; idx_cntr < (this->d_centro).size();
          idx_cntr++) {
        sutra_wfs *wfs = this->d_centro[idx_cntr]->wfs;
        float tmp_noise = wfs->noise;
        wfs->noise = -1;
        wfs->kernconv = true;
        wfs->sensor_trace(ydm, 1);
        wfs->comp_image();
        wfs->noise = tmp_noise;
        wfs->kernconv = false;
      }
      device = this->d_control[ncntrl]->d_centroids->getDevice();
      do_centroids(ncntrl, true);
      convert_centro(*this->d_control[ncntrl]->d_centroids,
          *this->d_control[ncntrl]->d_centroids, 0, 0.5f / p->second->push4imat,
          this->d_control[ncntrl]->d_centroids->getNbElem(),
          current_context->get_device(device));

      float alphai = -1.0f;
      cublasSaxpy(current_context->get_cublasHandle(),
          this->d_control[ncntrl]->d_centroids->getNbElem(), &alphai,
          this->d_control[ncntrl]->d_centroids->getData(), 1, (*d_imat)[inds1],
          1);

      p->second->reset_shape();
      inds1 += this->d_control[ncntrl]->nslope();
    }
    p++;
  }

  return EXIT_SUCCESS;
}

int sutra_rtc::do_imat_geom(int ncntrl, sutra_dms *ydm,
    int type) {
  current_context->set_activeDevice(device,1);
  map<type_screen, sutra_dm *>::iterator p;
  p = this->d_control[ncntrl]->d_dmseen.begin();
  int inds1, inds2;
  inds1 = 0;
  while (p != this->d_control[ncntrl]->d_dmseen.end()) {
    for (int j = 0; j < p->second->ninflu; j++) {
      p->second->comp_oneactu(j, p->second->push4imat); //
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
        if(wfs->type != "sh"){
          sutra_wfs_sh *_wfs = dynamic_cast<sutra_wfs_sh *>(wfs);
          _wfs->slopes_geom(type, d_imat[inds1 + inds2]);
        } else if(wfs->type != "geo"){
          sutra_wfs_geom *_wfs = dynamic_cast<sutra_wfs_geom *>(wfs);
          _wfs->slopes_geom(type, d_imat[inds1 + inds2]);
        } else {
          DEBUG_TRACE("wfs should be a SH or a geo");
          return EXIT_FAILURE;
        }
        inds2 += 2 * wfs->nvalid;
      }
      p->second->reset_shape();
      inds1 += this->d_control[ncntrl]->nslope();
    }
    p++;
  }
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
  return do_centroids(ncntrl, false);
}

int sutra_rtc::do_centroids(int ncntrl, bool imat) {
  current_context->set_activeDevice(device,1);
  int inds2 = 0;

  for (size_t idx_cntr = 0; idx_cntr < (this->d_centro).size(); idx_cntr++) {

    this->d_centro[idx_cntr]->get_cog((*this->d_control[ncntrl]->d_centroids)[inds2]);

    inds2 += 2 * this->d_centro[idx_cntr]->wfs->nvalid;
  }

  return EXIT_SUCCESS;
}

int sutra_rtc::do_centroids_geom(int ncntrl) {
  current_context->set_activeDevice(device,1);
  int inds2 = 0;

  for (size_t idx_cntr = 0; idx_cntr < (this->d_centro).size(); idx_cntr++) {

    sutra_wfs *wfs = this->d_centro[idx_cntr]->wfs;
    if(wfs->type == "sh"){
      sutra_wfs_sh *_wfs = dynamic_cast<sutra_wfs_sh *>(wfs);
      _wfs->slopes_geom(0,(*this->d_control[ncntrl]->d_centroids)[inds2]);

    } else if(wfs->type == "geo"){
      sutra_wfs_geom *_wfs = dynamic_cast<sutra_wfs_geom *>(wfs);
      _wfs->slopes_geom(0,(*this->d_control[ncntrl]->d_centroids)[inds2]);

    } else {
      DEBUG_TRACE("wfs could be a SH or a geo");
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

int sutra_rtc::do_control(int ncntrl) {
  current_context->set_activeDevice(device,1);

  if(this->d_control[ncntrl]->open_loop){
	  return EXIT_SUCCESS;
  }
  else
	  this->d_control[ncntrl]->comp_com();

  return EXIT_SUCCESS;
}

int sutra_rtc::apply_control(int ncntrl, sutra_dms *ydm) {
  current_context->set_activeDevice(device,1);

  this->d_control[ncntrl]->comp_voltage();

  map<type_screen, sutra_dm *>::iterator p;
  p = this->d_control[ncntrl]->d_dmseen.begin();
  int idx = 0;
  if ((this->d_control[ncntrl]->get_type().compare("ls") == 0)
      || (this->d_control[ncntrl]->get_type().compare("mv") == 0)
      || (this->d_control[ncntrl]->get_type().compare("geo") == 0)) {
    // "streamed" controllers case

    while (p != this->d_control[ncntrl]->d_dmseen.end()) {
      int nstreams = this->d_control[ncntrl]->streams->get_nbStreams();
      if (nstreams > p->second->ninflu) {
        for (int i = 0; i < nstreams; i++) {
          int istart = i * p->second->ninflu / nstreams;
          cutilSafeCall(
              cudaMemcpyAsync((*p->second->d_comm)[istart],
                  (*this->d_control[ncntrl]->d_voltage)[idx + istart],
                  sizeof(float) * p->second->ninflu / nstreams,
                  cudaMemcpyDeviceToDevice,
                  (*this->d_control[ncntrl]->streams)[i]));
          p->second->comp_shape();
        }
      } else {
        p->second->comp_shape((*this->d_control[ncntrl]->d_voltage)[idx]);
      }
      idx += p->second->ninflu;
      p++;
    }
  } else { // "non-streamed" controllers
    while (p != this->d_control[ncntrl]->d_dmseen.end()) {
      p->second->comp_shape((*this->d_control[ncntrl]->d_voltage)[idx]);
      idx += p->second->ninflu;
      p++;
    }
  }

  return EXIT_SUCCESS;
}

