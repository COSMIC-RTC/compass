/**
 * \file sutra_wfs_pyr_pyrhr.h
 *
 * \class sutra_wfs_pyr_pyrhr
 *
 * \ingroup libsutra
 *
 * \brief this class provides the wfs_pyr_pyrhr features to COMPASS
 *
 * \authors Damien Gratadour & Arnaud Sevin & Florian Ferreira
 *
 * \version 4.2.0
 *
 * \date 2011/01/28
 *
 */
#ifndef _SUTRA_WFS_PYR_PYRHR_H_
#define _SUTRA_WFS_PYR_PYRHR_H_

#include <sutra_lgs.h>
#include <sutra_phase.h>
#include <sutra_target.h>
#include <sutra_telemetry.h>
#include <sutra_telescope.h>
#include <sutra_wfs.h>
#include <map>
#include <vector>

using std::string;
class sutra_wfs_pyr_pyrhr : public sutra_wfs {
 public:
  long npupils;
  bool compute_pyrfocalplane;
  carma_obj<float> *d_hrimg;
  carma_obj<float> *d_submask;
  carma_obj<float> *d_psum;
  carma_obj<float> *d_pyrfocalplane;
  carma_obj<cuFloatComplex> *d_phalfxy;
  carma_obj<cuFloatComplex> *d_poffsets;

  carma_host_obj<float> *pyr_cx;
  carma_host_obj<float> *pyr_cy;
  carma_host_obj<float> *pyr_mod_weights;

 public:
  sutra_wfs_pyr_pyrhr(carma_context *context, sutra_telescope *d_tel,
                      carma_obj<cuFloatComplex> *d_camplipup,
                      carma_obj<cuFloatComplex> *d_camplifoc,
                      carma_obj<cuFloatComplex> *d_fttotim, long nxsub,
                      long nvalid, long npupils, long npix, long nphase,
                      long nrebin, long nfft, long ntot, long npup, float pdiam,
                      float nphotons, float nphot4imat, int lgs, bool fakecam,
                      int maxFluxPerPix, int maxPixValue, bool roket,
                      int device);
  sutra_wfs_pyr_pyrhr(carma_context *context, sutra_telescope *d_tel,
                      carma_obj<cuFloatComplex> *d_camplipup,
                      carma_obj<cuFloatComplex> *d_camplifoc,
                      carma_obj<cuFloatComplex> *d_fttotim, long nxsub,
                      long nvalid, long npupils, long npix, long nphase,
                      long nrebin, long nfft, long ntot, long npup, float pdiam,
                      float nphotons, float nphot4imat, int lgs, bool fakecam,
                      int maxFluxPerPix, int maxPixValue, bool roket,
                      int nbdevices, int *devices);
  ~sutra_wfs_pyr_pyrhr();

  int loadarrays(cuFloatComplex *halfxy, float *cx, float *cy, float *weights,
                 float *sincar, float *submask, int *validsubsx,
                 int *validsubsy, int *phasemap, float *fluxPerSub);
  int set_submask(float *submask);

  int fill_binimage(int async = 0);
  int comp_image(bool noise = true);
  void comp_modulation(int cpt);

  int copyValidPix(float *img, int *validx, int *validy, int im_dim);
  int set_pyr_modulation(float *cx, float *cy, int npts);
  int set_pyr_modulation(float *cx, float *cy, float *weights, int npts);
  int set_pyr_mod_weights(float *weights, int npts);

  int define_mpi_rank(int rank, int size) { return EXIT_SUCCESS; }
  int allocate_buffers(map<vector<int>, cufftHandle *> campli_plans,
                       map<vector<int>, cufftHandle *> fttotim_plans) {
    return EXIT_SUCCESS;
  }
  int comp_nphot(float ittime, float optthroughput, float diam, float cobs,
                 float zerop, float gsmag);

 private:
  int comp_generic();
  std::vector<carma_obj<cuFloatComplex> *> d_camplipup_ngpu;
  std::vector<carma_obj<cuFloatComplex> *> d_camplifoc_ngpu;
  std::vector<carma_obj<cuFloatComplex> *> d_phalfxy_ngpu;
  std::vector<carma_obj<cuFloatComplex> *> d_fttotim_ngpu;
  std::vector<carma_obj<float> *> d_pyrfocalplane_ngpu;
  std::vector<carma_obj<float> *> d_screen_ngpu;
  std::vector<carma_obj<float> *> d_hrimg_ngpu;
  std::vector<carma_obj<float> *> d_submask_ngpu;
};

#endif  // _SUTRA_WFS_PYR_PYRHR_H_
