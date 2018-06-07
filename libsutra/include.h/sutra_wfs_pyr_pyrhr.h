#ifndef _SUTRA_WFS_PYR_PYRHR_H_
#define _SUTRA_WFS_PYR_PYRHR_H_

#include <sutra_lgs.h>
#include <sutra_phase.h>
#include <sutra_target.h>
#include <sutra_telemetry.h>
#include <sutra_telescope.h>
#include <sutra_wfs_pyr.h>
#include <map>
#include <vector>

using std::string;
class sutra_sensors;
class sutra_wfs_pyr_pyrhr : public sutra_wfs_pyr {
 public:
 public:
  sutra_wfs_pyr_pyrhr(carma_context *context, sutra_telescope *d_tel,
                      sutra_sensors *sensors, long nxsub, long nvalid,
                      long npix, long nphase, long nrebin, long nfft, long ntot,
                      long npup, float pdiam, float nphotons, float nphot4imat,
                      int lgs, int device);
  sutra_wfs_pyr_pyrhr(carma_context *context, sutra_telescope *d_tel,
                      sutra_sensors *sensors, long nxsub, long nvalid,
                      long npix, long nphase, long nrebin, long nfft, long ntot,
                      long npup, float pdiam, float nphotons, float nphot4imat,
                      int lgs, int nbdevices, int *devices);
  sutra_wfs_pyr_pyrhr(const sutra_wfs_pyr_pyrhr &wfs);
  ~sutra_wfs_pyr_pyrhr();

  int wfs_initarrays(cuFloatComplex *halfxy, float *cx, float *cy,
                     float *sincar, float *submask, int *validsubsx,
                     int *validsubsy, int *phasemap, float *fluxPerSub);
  int set_submask(float *submask);

  int slopes_geom(int type, float *slopes);
  int slopes_geom(int type);

  int comp_image();
  void comp_modulation(int cpt);

  int copyValidPix(float *img, int *validx, int *validy, int im_dim);
  int set_pyr_modulation(float *cx, float *cy, int npts);

 private:
  int comp_generic();
  std::vector<carma_obj<cuFloatComplex> *> d_camplipup_ngpu;
  std::vector<carma_obj<cuFloatComplex> *> d_camplifoc_ngpu;
  std::vector<carma_obj<cuFloatComplex> *> d_phalfxy_ngpu;
  std::vector<carma_obj<cuFloatComplex> *> d_fttotim_ngpu;
  std::vector<carma_obj<float> *> d_screen_ngpu;
  std::vector<carma_obj<float> *> d_hrimg_ngpu;
  std::vector<carma_obj<float> *> d_submask_ngpu;
};

#endif  // _SUTRA_WFS_PYR_PYRHR_H_
