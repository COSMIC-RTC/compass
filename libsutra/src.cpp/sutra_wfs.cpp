#include <carma_utils.h>
#include <sutra_ao_utils.h>
#include <sutra_wfs.h>
#include <sutra_wfs_geom.h>
#include <sutra_wfs_pyr_pyr4.h>
#include <sutra_wfs_pyr_pyrhr.h>
#include <sutra_wfs_pyr_roof.h>
#include <sutra_wfs_sh.h>

int compute_nmaxhr(long nvalid) {
  // this is the big array => we use nmaxhr and treat it sequentially

  int mnmax = 500;
  int nmaxhr = mnmax;
  if (nvalid > 2 * mnmax) {
    int tmp0 = nvalid % mnmax;
    int tmp = 0;
    for (int cc = 1; cc < mnmax / 5; cc++) {
      tmp = nvalid % (mnmax + cc);
      if ((tmp > tmp0) || (tmp == 0)) {
        if (tmp == 0)
          tmp0 = 2 * mnmax;
        else
          tmp = tmp0;

        nmaxhr = mnmax + cc;
      }
    }
    return nmaxhr;
  }
  return nvalid;
}

sutra_wfs::sutra_wfs(carma_context *context, sutra_telescope *d_tel,
                     sutra_sensors *sensors, string type, long nxsub,
                     long nvalid, long npix, long nphase, long nrebin,
                     long nfft, long ntot, long npup, float pdiam,
                     float nphotons, float nphot4imat, int lgs, int device)
    : device(device),
      type(type),
      nxsub(nxsub),
      nvalid(nvalid),
      npix(npix),
      nrebin(nrebin),
      nfft(nfft),
      ntot(ntot),
      npup(npup),
      nphase(nphase),
      nmaxhr(nvalid),
      nffthr(1),
      subapd(pdiam),
      nphot(nphotons),
      nphot4imat(nphot4imat),
      noise(0),
      lgs(lgs),
      kernconv(false),
      error_budget(sensors->error_budget),
      campli_plan(nullptr),
      fttotim_plan(nullptr),
      d_ftkernel(nullptr),
      d_camplipup(nullptr),
      d_camplifoc(nullptr),
      d_fttotim(nullptr),
      d_pupil(d_tel->d_pupil_m),
      d_bincube(nullptr),
      d_bincube_notnoisy(nullptr),
      d_binimg(nullptr),
      d_binimg_notnoisy(nullptr),
      d_subsum(nullptr),
      d_offsets(nullptr),
      d_fluxPerSub(nullptr),
      d_sincar(nullptr),
      d_hrmap(nullptr),
      d_slopes(nullptr),
      image_telemetry(nullptr),
      d_gs(nullptr),
      streams(nullptr),
      nstreams(0),
      d_phasemap(nullptr),
      d_validsubsx(nullptr),
      d_validsubsy(nullptr),
      current_context(context),
      offset(0),
      nvalid_tot(nvalid),
      rank(0),
      displ_bincube(nullptr),
      count_bincube(nullptr) {
  if (sensors != nullptr) {
    this->d_camplipup = sensors->d_camplipup;
    this->d_camplifoc = sensors->d_camplifoc;
    this->d_fttotim = sensors->d_fttotim;
  }
}

int sutra_wfs::wfs_initgs(sutra_sensors *sensors, float xpos, float ypos,
                          float lambda, float mag, float zerop, long size,
                          float noise, long seed, float G, float thetaML,
                          float dx, float dy) {
  current_context->set_activeDevice(device, 1);
  this->d_gs = new sutra_source(current_context, xpos, ypos, lambda, mag, zerop,
                                size, "wfs", this->device);
  this->d_gs->G = G;
  this->d_gs->thetaML = thetaML;
  this->d_gs->dx = dx;
  this->d_gs->dy = dy;

  set_noise(noise, seed);

  if (this->lgs) {
    this->d_gs->d_lgs = new sutra_lgs(current_context, sensors, this->nvalid,
                                      this->ntot, this->nmaxhr);
    this->d_gs->lgs = this->lgs;
  }

  return EXIT_SUCCESS;
}

int sutra_wfs::set_noise(float noise, long seed) {
  current_context->set_activeDevice(device, 1);
  this->noise = noise;
  if (this->type != "pyrhr") {
    if (noise > -1) {
      this->d_bincube->init_prng(seed);
      // this->d_bincube->prng('N', noise, 0.0f);
    }
    if (noise > 0) {
      this->d_binimg->init_prng(seed);
      // this->d_binimg->prng('N', noise, 0.0f);
    }
  } else {
    if (noise > -1) {
      this->d_binimg->init_prng(seed);
    }
  }
  return EXIT_SUCCESS;
}

int sutra_wfs::load_pupil(float *pupil) {
  current_context->set_activeDevice(device, 1);
  this->d_pupil->host2device(pupil);
  int nbdevices = d_pupil_ngpu.size();
  for (int ndevice = 1; ndevice < nbdevices; ndevice++) {
    current_context->set_activeDevice(ndevice, 1);
    d_pupil_ngpu[ndevice]->host2device(pupil);
  }
  current_context->set_activeDevice(device, 1);
  return EXIT_SUCCESS;
}

int sutra_wfs::load_kernels(float *lgskern) {
  current_context->set_activeDevice(device, 1);
  if (this->lgs)
    this->d_gs->d_lgs->load_kernels(lgskern,
                                    this->current_context->get_device(device));

  return EXIT_SUCCESS;
}

int sutra_wfs::sensor_trace(sutra_atmos *yatmos) {
  current_context->set_activeDevice(device, 1);
  // do raytracing to get the phase
  this->d_gs->raytrace(yatmos);

  // same with textures
  // apparently not working for large sizes...
  // this->d_gs->raytrace_shm(yatmos);

  return EXIT_SUCCESS;
}

int sutra_wfs::sensor_trace(sutra_dms *ydm, int rst) {
  current_context->set_activeDevice(device, 1);
  // do raytracing to get the phase
  this->d_gs->raytrace(ydm, rst, 0);

  return EXIT_SUCCESS;
}

int sutra_wfs::sensor_trace(int rst) {
  current_context->set_activeDevice(device, 1);
  // do raytracing to get the phase
  this->d_gs->raytrace(rst);

  return EXIT_SUCCESS;
}

int sutra_wfs::sensor_trace(sutra_atmos *yatmos, sutra_dms *ydms) {
  this->d_gs->raytrace(yatmos);
  this->d_gs->raytrace(ydms, 0, 0);

  return EXIT_SUCCESS;
}

sutra_sensors::sutra_sensors(carma_context *context, sutra_telescope *d_tel,
                             char **type, int nwfs, long *nxsub, long *nvalid,
                             long *npix, long *nphase, long *nrebin, long *nfft,
                             long *ntot, long *npup, float *pdiam, float *nphot,
                             float *nphot4imat, int *lgs, int device,
                             bool error_budget)
    : device(device),
      error_budget(error_budget),
      current_context(context),
      d_camplipup(nullptr),
      d_camplifoc(nullptr),
      d_fttotim(nullptr),
      d_ftlgskern(nullptr),
      d_lgskern(nullptr) {
  current_context->set_activeDevice(device, 1);
  // DEBUG_TRACE("Before create sensors : ");printMemInfo();
  if (strcmp(type[0], "sh") == 0) {
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
    if (strcmp(type[i], "sh") == 0)
      wfs = new sutra_wfs_sh(context, d_tel, this, nxsub[i], nvalid[i], npix[i],
                             nphase[i], nrebin[i], nfft[i], ntot[i], npup[i],
                             pdiam[i], nphot[i], nphot4imat[i], lgs[i], device);
    if (strcmp(type[i], "pyr") == 0)
      wfs = new sutra_wfs_pyr_pyr4(context, d_tel, this, nxsub[i], nvalid[i],
                                   npix[i], nphase[i], nrebin[i], nfft[i],
                                   ntot[i], npup[i], pdiam[i], nphot[i],
                                   nphot4imat[i], lgs[i], device);
    if (strcmp(type[i], "pyrhr") == 0) {
      const int ngpu = context->get_ndevice();
      DEBUG_TRACE("using pyrhr with %d GPUs", ngpu);
      if (ngpu == 1) {
        wfs = new sutra_wfs_pyr_pyrhr(context, d_tel, this, nxsub[i], nvalid[i],
                                      npix[i], nphase[i], nrebin[i], nfft[i],
                                      ntot[i], npup[i], pdiam[i], nphot[i],
                                      nphot4imat[i], lgs[i], device);
      } else {
        int devices[ngpu];
        for (int i = 0; i < ngpu; i++) {
          devices[i] = i;
        }
        wfs = new sutra_wfs_pyr_pyrhr(context, d_tel, this, nxsub[i], nvalid[i],
                                      npix[i], nphase[i], nrebin[i], nfft[i],
                                      ntot[i], npup[i], pdiam[i], nphot[i],
                                      nphot4imat[i], lgs[i], ngpu, devices);
      }
    }
    if (strcmp(type[i], "roof") == 0)
      wfs = new sutra_wfs_pyr_roof(context, d_tel, this, nxsub[i], nvalid[i],
                                   npix[i], nphase[i], nrebin[i], nfft[i],
                                   ntot[i], npup[i], pdiam[i], nphot[i],
                                   nphot4imat[i], lgs[i], device);
    d_wfs.push_back(wfs);

    //		DEBUG_TRACE("After creating wfs #%d : ",i);printMemInfo();
  }
  //	DEBUG_TRACE("Final sensors : ");printMemInfo();
}

sutra_sensors::sutra_sensors(carma_context *context, sutra_telescope *d_tel,
                             int nwfs, long *nxsub, long *nvalid, long *nphase,
                             long npup, float *pdiam, int device)
    : device(device),
      error_budget(false),
      current_context(context),
      d_camplipup(nullptr),
      d_camplifoc(nullptr),
      d_fttotim(nullptr),
      d_ftlgskern(nullptr),
      d_lgskern(nullptr) {
  this->current_context = context;
  this->device = device;
  current_context->set_activeDevice(device, 1);
  DEBUG_TRACE("device %d", device);

  for (int i = 0; i < nwfs; i++) {
    d_wfs.push_back(new sutra_wfs_geom(context, d_tel, nxsub[i], nvalid[i],
                                       nphase[i], npup, pdiam[i], device));
  }
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
    this->d_wfs[idx]->allocate_buffers(this);
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

int sutra_sensors::sensors_initgs(float *xpos, float *ypos, float *lambda,
                                  float *mag, float zerop, long *size,
                                  float *noise, long *seed, float *G,
                                  float *thetaML, float *dx, float *dy) {
  current_context->set_activeDevice(device, 1);
  for (size_t idx = 0; idx < this->d_wfs.size(); idx++) {
    this->d_wfs[idx]->wfs_initgs(
        this, xpos[idx], ypos[idx], lambda[idx], mag[idx], zerop, size[idx],
        noise[idx], seed[idx], G[idx], thetaML[idx], dx[idx], dy[idx]);
  }
  return EXIT_SUCCESS;
}
int sutra_sensors::sensors_initgs(float *xpos, float *ypos, float *lambda,
                                  float *mag, float zerop, long *size,
                                  float *noise, float *G, float *thetaML,
                                  float *dx, float *dy) {
  current_context->set_activeDevice(device, 1);
  for (size_t idx = 0; idx < this->d_wfs.size(); idx++) {
    this->d_wfs[idx]->wfs_initgs(
        this, xpos[idx], ypos[idx], lambda[idx], mag[idx], zerop, size[idx],
        noise[idx], 1234 + idx, G[idx], thetaML[idx], dx[idx], dy[idx]);
  }
  return EXIT_SUCCESS;
}
int sutra_sensors::sensors_initgs(float *xpos, float *ypos, float *lambda,
                                  float *mag, float zerop, long *size, float *G,
                                  float *thetaML, float *dx, float *dy) {
  current_context->set_activeDevice(device, 1);
  for (size_t idx = 0; idx < this->d_wfs.size(); idx++) {
    this->d_wfs[idx]->wfs_initgs(this, xpos[idx], ypos[idx], lambda[idx],
                                 mag[idx], zerop, size[idx], -1, 1234 + idx,
                                 G[idx], thetaML[idx], dx[idx], dy[idx]);
  }
  return EXIT_SUCCESS;
}

int sutra_sensors::set_noise(int nwfs, float noise, long seed) {
  this->d_wfs[nwfs]->set_noise(noise, seed);
  return EXIT_SUCCESS;
}
