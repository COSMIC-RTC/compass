#include <sutra_sensors.h>

sutra_sensors::sutra_sensors(carma_context *context, sutra_telescope *d_tel,
                             vector<string> type, int nwfs, long *nxsub,
                             long *nvalid, long *npix, long *nphase,
                             long *nrebin, long *nfft, long *ntot, long *npup,
                             float *pdiam, float *nphot, float *nphot4imat,
                             int *lgs, int device, bool roket)
    : device(device),
      current_context(context),
      d_camplipup(nullptr),
      d_camplifoc(nullptr),
      d_fttotim(nullptr),
      d_ftlgskern(nullptr),
      d_lgskern(nullptr) {
  current_context->set_activeDevice(device, 1);
  // DEBUG_TRACE("Before create sensors : ");printMemInfo();
  if (type[0].compare("sh") == 0) {
    int maxnfft = nfft[0];
    int maxntot = ntot[0];
    int maxnvalid = 0;
    int maxnvalid_tot = nvalid[0];
    int is_lgs = (lgs[0] > 0 ? 1 : 0);
    for (int i = 1; i < nwfs; i++) {
      if (nvalid[i] > maxnvalid_tot) maxnvalid_tot = nvalid[i];
      if (ntot[i] > maxntot) {
        maxntot = ntot[i];
      }
      if (nfft[i] > maxnfft) {
        maxnfft = nfft[i];
      }
      if (ntot[i] == nfft[i]) {
        if (nvalid[i] > maxnvalid) {
          maxnvalid = nvalid[i];
        }
      }
      if (lgs[i] > 0) is_lgs = 1;
    }
    // DEBUG_TRACE("maxntot : %d maxnfft : %d maxnvalid : %d nmaxhr : %d \n
    // ",maxntot,maxnfft,maxnvalid,compute_nmaxhr(nvalid[wfs4ntot]));
    long dims_data3[4] = {
        3, maxnfft, maxnfft, maxnvalid_tot /*nvalid[wfs4nfft]*/
    };
    this->d_camplifoc = new carma_obj<cuFloatComplex>(context, dims_data3);
    this->d_camplipup = new carma_obj<cuFloatComplex>(context, dims_data3);

    dims_data3[1] = maxntot;
    dims_data3[2] = maxntot;
    dims_data3[3] = compute_nmaxhr(maxnvalid_tot /*nvalid[wfs4ntot]*/);
    if (maxnvalid > dims_data3[3]) dims_data3[3] = maxnvalid;
    this->d_fttotim = new carma_obj<cuFloatComplex>(context, dims_data3);
    if (is_lgs) {
      this->d_ftlgskern = new carma_obj<cuFloatComplex>(context, dims_data3);
      this->d_lgskern = new carma_obj<float>(context, dims_data3);
    } else {
      this->d_ftlgskern = 0L;
      this->d_lgskern = 0L;
    }
  } else {  // Pyramid case : shared arrays are not used
    this->d_camplifoc = 0L;
    this->d_camplipup = 0L;
    this->d_fttotim = 0L;
    this->d_ftlgskern = 0L;
    this->d_lgskern = 0L;
  }

  // DEBUG_TRACE("After creating sensors arrays : ");printMemInfo();
  for (int i = 0; i < nwfs; i++) {
    sutra_wfs *wfs = NULL;
    if (type[i].compare("sh") == 0)
      wfs = new sutra_wfs_sh(
          context, d_tel, this->d_camplipup, this->d_camplifoc, this->d_fttotim,
          nxsub[i], nvalid[i], npix[i], nphase[i], nrebin[i], nfft[i], ntot[i],
          npup[i], pdiam[i], nphot[i], nphot4imat[i], lgs[i], roket, device);

    if (type[i].compare("pyrhr") == 0) {
      const int ngpu = context->get_ndevice();
      DEBUG_TRACE("using pyrhr with %d GPUs", ngpu);
      if (ngpu == 1) {
        wfs = new sutra_wfs_pyr_pyrhr(
            context, d_tel, this->d_camplipup, this->d_camplifoc,
            this->d_fttotim, nxsub[i], nvalid[i], npix[i], nphase[i], nrebin[i],
            nfft[i], ntot[i], npup[i], pdiam[i], nphot[i], nphot4imat[i],
            lgs[i], roket, device);
      } else {
        int devices[ngpu];
        for (int i = 0; i < ngpu; i++) {
          devices[i] = i;
        }
        wfs = new sutra_wfs_pyr_pyrhr(
            context, d_tel, this->d_camplipup, this->d_camplifoc,
            this->d_fttotim, nxsub[i], nvalid[i], npix[i], nphase[i], nrebin[i],
            nfft[i], ntot[i], npup[i], pdiam[i], nphot[i], nphot4imat[i],
            lgs[i], ngpu, roket, devices);
      }
    }

    d_wfs.push_back(wfs);

    //		DEBUG_TRACE("After creating wfs #%d : ",i);printMemInfo();
  }
  //	DEBUG_TRACE("Final sensors : ");printMemInfo();
  this->allocate_buffers();
}

sutra_sensors::~sutra_sensors() {
  current_context->set_activeDevice(device, 1);
  //  for (size_t idx = 0; idx < (this->d_wfs).size(); idx++) {
  while ((this->d_wfs).size() > 0) {
    delete this->d_wfs.back();
    d_wfs.pop_back();
  }
  map<vector<int>, cufftHandle *>::iterator it;
  for (it = campli_plans.begin(); it != campli_plans.end(); it++) {
    cufftDestroy(*it->second);
    free(it->second);
  }
  for (it = fttotim_plans.begin(); it != fttotim_plans.end(); it++) {
    cufftDestroy(*it->second);
    free(it->second);
  }
  for (it = ftlgskern_plans.begin(); it != ftlgskern_plans.end(); it++) {
    cufftDestroy(*it->second);
    free(it->second);
  }

  if (this->d_camplifoc != nullptr) delete this->d_camplifoc;
  if (this->d_camplipup != nullptr) delete this->d_camplipup;
  if (this->d_fttotim != nullptr) delete this->d_fttotim;
  if (this->d_lgskern != nullptr) delete this->d_lgskern;
  if (this->d_ftlgskern != nullptr) delete this->d_ftlgskern;
}

int sutra_sensors::allocate_buffers() {
  current_context->set_activeDevice(device, 1);
  for (size_t idx = 0; idx < this->d_wfs.size(); idx++) {
    this->d_wfs[idx]->allocate_buffers(this->campli_plans, this->fttotim_plans);
  }
  return EXIT_SUCCESS;
}

int sutra_sensors::define_mpi_rank(int rank, int size) {
  current_context->set_activeDevice(device, 1);
  for (size_t idx = 0; idx < this->d_wfs.size(); idx++) {
    this->d_wfs[idx]->define_mpi_rank(rank, size);
  }
  return EXIT_SUCCESS;
}

int sutra_sensors::initgs(float *xpos, float *ypos, float *lambda, float *mag,
                          float zerop, long *size, float *noise, long *seed,
                          float *G, float *thetaML, float *dx, float *dy) {
  current_context->set_activeDevice(device, 1);
  for (size_t idx = 0; idx < this->d_wfs.size(); idx++) {
    this->d_wfs[idx]->wfs_initgs(
        this->d_lgskern, this->d_ftlgskern, this->ftlgskern_plans, xpos[idx],
        ypos[idx], lambda[idx], mag[idx], zerop, size[idx], noise[idx],
        seed[idx], G[idx], thetaML[idx], dx[idx], dy[idx]);
  }
  return EXIT_SUCCESS;
}
int sutra_sensors::initgs(float *xpos, float *ypos, float *lambda, float *mag,
                          float zerop, long *size, float *noise, float *G,
                          float *thetaML, float *dx, float *dy) {
  current_context->set_activeDevice(device, 1);
  for (size_t idx = 0; idx < this->d_wfs.size(); idx++) {
    this->d_wfs[idx]->wfs_initgs(
        this->d_lgskern, this->d_ftlgskern, this->ftlgskern_plans, xpos[idx],
        ypos[idx], lambda[idx], mag[idx], zerop, size[idx], noise[idx],
        1234 + idx, G[idx], thetaML[idx], dx[idx], dy[idx]);
  }
  return EXIT_SUCCESS;
}
int sutra_sensors::initgs(float *xpos, float *ypos, float *lambda, float *mag,
                          float zerop, long *size, float *G, float *thetaML,
                          float *dx, float *dy) {
  current_context->set_activeDevice(device, 1);
  for (size_t idx = 0; idx < this->d_wfs.size(); idx++) {
    this->d_wfs[idx]->wfs_initgs(
        this->d_lgskern, this->d_ftlgskern, this->ftlgskern_plans, xpos[idx],
        ypos[idx], lambda[idx], mag[idx], zerop, size[idx], -1, 1234 + idx,
        G[idx], thetaML[idx], dx[idx], dy[idx]);
  }
  return EXIT_SUCCESS;
}

int sutra_sensors::set_noise(int nwfs, float noise, long seed) {
  this->d_wfs[nwfs]->set_noise(noise, seed);
  return EXIT_SUCCESS;
}
