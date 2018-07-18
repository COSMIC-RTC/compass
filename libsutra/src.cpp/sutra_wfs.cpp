#include <sutra_wfs.h>

sutra_wfs::sutra_wfs(carma_context *context, sutra_telescope *d_tel,
                     carma_obj<cuFloatComplex> *d_camplipup,
                     carma_obj<cuFloatComplex> *d_camplifoc,
                     carma_obj<cuFloatComplex> *d_fttotim, string type,
                     long nxsub, long nvalid, long npix, long nphase,
                     long nrebin, long nfft, long ntot, long npup, float pdiam,
                     float nphotons, float nphot4imat, int lgs, bool roket,
                     int device)
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
      roket(roket),
      campli_plan(nullptr),
      fttotim_plan(nullptr),
      d_ftkernel(nullptr),
      d_camplipup(d_camplipup),
      d_camplifoc(d_camplifoc),
      d_fttotim(d_fttotim),
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
      count_bincube(nullptr) {}

int sutra_wfs::wfs_initgs(carma_obj<float> *d_lgskern,
                          carma_obj<cuFloatComplex> *d_ftlgskern,
                          map<vector<int>, cufftHandle *> ftlgskern_plans,
                          float xpos, float ypos, float lambda, float mag,
                          float zerop, long size, float noise, long seed,
                          float G, float thetaML, float dx, float dy) {
  current_context->set_activeDevice(device, 1);
  this->d_gs = new sutra_source(current_context, xpos, ypos, lambda, mag, zerop,
                                size, "wfs", this->device);
  this->d_gs->G = G;
  this->d_gs->thetaML = thetaML;
  this->d_gs->dx = dx;
  this->d_gs->dy = dy;

  set_noise(noise, seed);

  if (this->lgs) {
    this->d_gs->d_lgs =
        new sutra_lgs(current_context, d_lgskern, d_ftlgskern, ftlgskern_plans,
                      this->nvalid, this->ntot, this->nmaxhr);
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

int sutra_wfs::set_pupil(float *pupil) {
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
