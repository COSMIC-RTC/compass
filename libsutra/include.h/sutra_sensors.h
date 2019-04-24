/**
 * \file sutra_sensors.h
 *
 * \class sutra_sensors
 *
 * \ingroup libsutra
 *
 * \brief this class provides the sensors features to COMPASS
 *
 * \authors Damien Gratadour & Arnaud Sevin & Florian Ferreira
 *
 * \version 1.0
 *
 * \date 2011/01/28
 *
 */
#ifndef _SUTRA_SENSORS_H_
#define _SUTRA_SENSORS_H_

#include <carma_utils.h>
#include <sutra_telescope.h>
#include <sutra_utils.h>
#include <sutra_wfs.h>
#include <sutra_wfs_pyr_pyrhr.h>
#include <sutra_wfs_sh.h>

#include <map>
#include <vector>

using std::map;
using std::string;
using std::vector;

class sutra_sensors {
 public:
  int device;
  bool roket;
  carma_context *current_context;
  size_t nsensors() { return d_wfs.size(); }
  vector<sutra_wfs *> d_wfs;
  map<vector<int>, cufftHandle *> campli_plans;
  map<vector<int>, cufftHandle *> fttotim_plans;
  map<vector<int>, cufftHandle *> ftlgskern_plans;

  carma_obj<cuFloatComplex> *d_camplipup;
  carma_obj<cuFloatComplex> *d_camplifoc;
  carma_obj<cuFloatComplex> *d_fttotim;
  carma_obj<cuFloatComplex> *d_ftlgskern;
  carma_obj<float> *d_lgskern;

 public:
  sutra_sensors(carma_context *context, sutra_telescope *d_tel,
                vector<string> type, int nwfs, long *nxsub, long *nvalid,
                long *npupils, long *npix, long *nphase, long *nrebin,
                long *nfft, long *ntot, long *npup, float *pdiam, float *nphot,
                float *nphot4imat, int *lgs, bool *fakecam, int *maxFluxPerPix,
                int *maxPixValue, int device, bool roket);
  ~sutra_sensors();

  int allocate_buffers();
  int define_mpi_rank(int rank, int size);
  int set_noise(int nwfs, float noise, long seed);

  int initgs(float *xpos, float *ypos, float *lambda, float *mag, float zerop,
             long *size, float *noise, long *seed, float *G, float *thetaML,
             float *dx, float *dy);
  int initgs(float *xpos, float *ypos, float *lambda, float *mag, float zerop,
             long *size, float *noise, float *G, float *thetaML, float *dx,
             float *dy);
  int initgs(float *xpos, float *ypos, float *lambda, float *mag, float zerop,
             long *size, float *G, float *thetaML, float *dx, float *dy);
};

// General utilities

#endif  // _SUTRA_SENSORS_H_
