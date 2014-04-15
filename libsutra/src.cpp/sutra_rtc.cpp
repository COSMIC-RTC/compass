#include <sutra_rtc.h>

sutra_rtc::sutra_rtc(carma_context *context) {
  this->current_context = context;
  this->device = context->get_activeDevice();
  this->delay = 0;
}

sutra_rtc::~sutra_rtc() {
  for (size_t idx = 0; idx < (this->d_centro).size(); idx++) {
    delete this->d_centro[(this->d_centro).size() - 1];
    d_centro.pop_back();
    ;
  }
  for (size_t idx = 0; idx < (this->d_control).size(); idx++) {
    delete this->d_control[(this->d_control).size() - 1];
    d_control.pop_back();
  }

  //delete this->current_context;
}

int sutra_rtc::add_centroider(long nwfs, long nvalid, float offset, float scale,
    long device, char *typec) {
  if (strcmp(typec, "bpcog") == 0)
    d_centro.push_back(
        new sutra_centroider_bpcog(current_context, nwfs, nvalid, offset, scale,
            device, 10));
  else if (strcmp(typec, "cog") == 0)
    d_centro.push_back(
        new sutra_centroider_cog(current_context, nwfs, nvalid, offset, scale,
            device));
  else if (strcmp(typec, "corr") == 0)
    d_centro.push_back(
        new sutra_centroider_corr(current_context, nwfs, nvalid, offset, scale,
            device));
  else if (strcmp(typec, "pyr") == 0)
    d_centro.push_back(
        new sutra_centroider_pyr(current_context, nwfs, nvalid, offset, scale,
            device));
  else if (strcmp(typec, "tcog") == 0)
    d_centro.push_back(
        new sutra_centroider_tcog(current_context, nwfs, nvalid, offset, scale,
            device));
  else if (strcmp(typec, "wcog") == 0)
    d_centro.push_back(
        new sutra_centroider_wcog(current_context, nwfs, nvalid, offset, scale,
            device));
  else {
    cerr << "centroider unknown\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

int sutra_rtc::add_controller(long nactu, long delay, long device,
    const char *typec) {
  int ncentroids = 0;
  for (size_t idx = 0; idx < (this->d_centro).size(); idx++)
    ncentroids += this->d_centro[idx]->nvalid;

  current_context->set_activeDevice(device);
  string type_ctr(typec);
  if (type_ctr.compare("ls") == 0) {
    d_control.push_back(
        new sutra_controller_ls(current_context, ncentroids, nactu, delay));
  } else if (type_ctr.compare("cured") == 0) {
    d_control.push_back(
        new sutra_controller_cured(current_context, ncentroids, nactu, delay));
  } else if (type_ctr.compare("mv") == 0) {
    d_control.push_back(
        new sutra_controller_mv(current_context, ncentroids, nactu, delay));
  } else if (type_ctr.compare("kalman") == 0) {
    carma_obj<float>* cD_Mo = calculate_D_Mo(current_context, ncentroids * 2);
    carma_obj<float>* cN_Act = calculate_N_Act(current_context, nactu);
    carma_obj<float>* cPROJ = calculate_btur();
    d_control.push_back(
        new sutra_controller_kalman(current_context, cD_Mo, cN_Act, cPROJ,
            false));
    if (cD_Mo)
      delete cD_Mo;
    if (cN_Act)
      delete cN_Act;
    if (cPROJ)
      delete cPROJ;

  } else {
    DEBUG_TRACE("Controller '%s' unknown\n", typec);
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

int sutra_rtc::rm_controller() {

  delete this->d_control[(this->d_control).size() - 1];
  d_control.pop_back();

  return EXIT_SUCCESS;
}

int sutra_rtc::do_imat(int ncntrl, sutra_sensors *sensors, sutra_dms *ydm) {
  carma_obj<float> *d_imat;
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
  p = ydm->d_dms.begin();
  int inds1;
  inds1 = 0;
  cout << "starting imat measurement" << endl;
  while (p != ydm->d_dms.end()) {
    for (int j = 0; j < p->second->ninflu; j++) {
      p->second->comp_oneactu(j, p->second->push4imat);
      for (size_t idx_cntr = 0; idx_cntr < (this->d_centro).size();
          idx_cntr++) {
        int nwfs = this->d_centro[idx_cntr]->nwfs;
        float tmp_noise = sensors->d_wfs[nwfs]->noise;
        sensors->d_wfs[nwfs]->noise = -1;
        sensors->d_wfs[nwfs]->kernconv = true;
        sensors->d_wfs[nwfs]->sensor_trace(ydm, 1);
        sensors->d_wfs[nwfs]->comp_image();
        sensors->d_wfs[nwfs]->noise = tmp_noise;
        sensors->d_wfs[nwfs]->kernconv = false;
      }
      //cout << "actu # " << j << endl;
      do_centroids(ncntrl, sensors, true);

      convert_centro(*this->d_control[ncntrl]->d_centroids,
          *this->d_control[ncntrl]->d_centroids, 0, 1.0f / p->second->push4imat,
          this->d_control[ncntrl]->d_centroids->getNbElem(),
          this->d_control[ncntrl]->d_centroids->getDevice());

      this->d_control[ncntrl]->d_centroids->copyInto((*d_imat)[inds1],
          this->d_control[ncntrl]->nslope());

      p->second->reset_shape();
      inds1 += this->d_control[ncntrl]->nslope();
    }
    p++;
  }

  return EXIT_SUCCESS;
}

int sutra_rtc::do_imat_geom(int ncntrl, sutra_sensors *sensors, sutra_dms *ydm,
    int type) {
  if (this->d_control[ncntrl]->get_type().compare("ls") == 0) {
    SCAST(sutra_controller_ls *, control, this->d_control[ncntrl]);
    map<type_screen, sutra_dm *>::iterator p;
    p = ydm->d_dms.begin();
    int inds1, inds2;
    inds1 = 0;
    while (p != ydm->d_dms.end()) {
      for (int j = 0; j < p->second->ninflu; j++) {
        p->second->comp_oneactu(j, p->second->push4imat); //
        inds2 = 0;
        for (size_t idx_cntr = 0; idx_cntr < (this->d_control).size();
            idx_cntr++) {
          int nwfs = this->d_centro[idx_cntr]->nwfs;

          sensors->d_wfs[nwfs]->sensor_trace(ydm, 1);
          //sensors->d_wfs[nwfs]->comp_image();

          sensors->d_wfs[nwfs]->slopes_geom(type,
              (*control->d_imat)[inds1 + inds2]);
          inds2 += 2 * sensors->d_wfs[nwfs]->nvalid;
        }
        p->second->reset_shape();
        inds1 += control->nslope();
      }
      p++;
    }
  }
  return EXIT_SUCCESS;
}

int sutra_rtc::do_centroids(sutra_sensors *sensors) {
  for (size_t idx_cntr = 0; idx_cntr < (this->d_centro).size(); idx_cntr++) {
    int nwfs = this->d_centro[idx_cntr]->nwfs;

    this->d_centro[idx_cntr]->get_cog(sensors->d_wfs[nwfs]);

  }

  return EXIT_SUCCESS;
}

int sutra_rtc::do_centroids(int ncntrl, sutra_sensors *sensors) {
  return do_centroids(ncntrl, sensors, false);
}

int sutra_rtc::do_centroids(int ncntrl, sutra_sensors *sensors, bool imat) {
  int inds2 = 0;

  for (size_t idx_cntr = 0; idx_cntr < (this->d_centro).size(); idx_cntr++) {

    int nwfs = this->d_centro[idx_cntr]->nwfs;

    this->d_centro[idx_cntr]->get_cog(sensors->d_wfs[nwfs],
        (*this->d_control[ncntrl]->d_centroids)[inds2]);

    inds2 += 2 * sensors->d_wfs[nwfs]->nvalid;
  }

  return EXIT_SUCCESS;
}

int sutra_rtc::do_control(int ncntrl, sutra_dms *ydm) {
  if (this->d_control[ncntrl]->get_type().compare("ls") == 0) {
    SCAST(sutra_controller_ls *, control, this->d_control[ncntrl]);
    //   fprintf(stderr, "[%s@%d] here!\n", __FILE__, __LINE__);

    control->frame_delay();
    control->comp_com();

    map<type_screen, sutra_dm *>::iterator p;
    p = ydm->d_dms.begin();
    int idx = 0;
    while (p != ydm->d_dms.end()) {
      int nstreams = control->streams->get_nbStreams();
      if (nstreams > p->second->ninflu) {
        for (int i = 0; i < nstreams; i++) {
          int istart = i * p->second->ninflu / nstreams;
          cutilSafeCall(
              cudaMemcpyAsync((*p->second->d_comm)[istart],
                  (*control->d_com)[idx + istart],
                  sizeof(float) * p->second->ninflu / nstreams,
                  cudaMemcpyDeviceToDevice, (*control->streams)[i]));
          p->second->comp_shape();
        }
      } else {
        p->second->comp_shape((*control->d_com)[idx]);
      }
      idx += p->second->ninflu;
      p++;
    }
  } else if (this->d_control[ncntrl]->get_type().compare("mv") == 0) {
    SCAST(sutra_controller_mv *, control, this->d_control[ncntrl]);
    //   fprintf(stderr, "[%s@%d] here!\n", __FILE__, __LINE__);

    control->frame_delay();
    control->comp_com();

    map<type_screen, sutra_dm *>::iterator p;
    p = ydm->d_dms.begin();
    int idx = 0;
    while (p != ydm->d_dms.end()) {
      int nstreams = control->streams->get_nbStreams();
      if (nstreams > p->second->ninflu) {
        for (int i = 0; i < nstreams; i++) {
          int istart = i * p->second->ninflu / nstreams;
          cutilSafeCall(
              cudaMemcpyAsync((*p->second->d_comm)[istart],
                  (*control->d_com)[idx + istart],
                  sizeof(float) * p->second->ninflu / nstreams,
                  cudaMemcpyDeviceToDevice, (*control->streams)[i]));
          p->second->comp_shape();
        }
      } else {
        /*
         float values[control->d_com->getNbElem()];
         control->d_com->device2host(values);
         for (int i=0 ; i<control->d_com->getNbElem() ; i++)
         cout << values[i] << " ";
         cout << endl;
         */
        p->second->comp_shape((*control->d_com)[idx]);
      }
      idx += p->second->ninflu;
      p++;
    }
  } else {
    this->d_control[ncntrl]->comp_com();
    map<type_screen, sutra_dm *>::iterator p;
    p = ydm->d_dms.begin();
    int idx = 0;
    while (p != ydm->d_dms.end()) {
      p->second->comp_shape((*this->d_control[ncntrl]->d_com)[idx]);
     p++;
    }
  }

  return EXIT_SUCCESS;
}

