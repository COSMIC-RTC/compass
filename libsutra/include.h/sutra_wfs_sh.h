#ifndef _SUTRA_WFS_SH_H_
#define _SUTRA_WFS_SH_H_

#include <sutra_lgs.h>
#include <sutra_phase.h>
#include <sutra_target.h>
#include <sutra_telemetry.h>
#include <sutra_wfs.h>
#include <map>
#include <vector>

class sutra_wfs_sh : public sutra_wfs {
 public:
  // sh only
  carma_obj<int> *d_binmap;
  carma_obj<int> *d_istart;  // nxsub
  carma_obj<int> *d_jstart;  // nxsub

 public:
  sutra_wfs_sh(carma_context *context, sutra_telescope *d_tel,
               sutra_sensors *sensors, long nxsub, long nvalid, long npix,
               long nphase, long nrebin, long nfft, long ntot, long npup,
               float pdiam, float nphotons, float nphot4imat, int lgs,
               int device);
  sutra_wfs_sh(const sutra_wfs_sh &wfs);
  ~sutra_wfs_sh();

  int define_mpi_rank(int rank, int size);
  int allocate_buffers(sutra_sensors *sensors);

  int wfs_initarrays(int *phasemap, int *hrmap, int *binmap, float *offsets,
                     float *fluxPerSub, int *validsubsx, int *validsubsy,
                     int *istart, int *jstart, cuFloatComplex *kernel);

  int slopes_geom(int type, float *slopes);
  int slopes_geom(int type);
  int fill_binimage(int async);
  int comp_image();

 private:
  int comp_generic();
};

#endif  // _SUTRA_WFS_SH_H_
