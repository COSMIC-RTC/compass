#ifndef _SUTRA_WFS_SH_H_
#define _SUTRA_WFS_SH_H_

#include <sutra_lgs.h>
#include <sutra_phase.h>
#include <sutra_target.h>
#include <sutra_telemetry.h>
#include <sutra_utils.h>
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
               carma_obj<cuFloatComplex> *d_camplipup,
               carma_obj<cuFloatComplex> *d_camplifoc,
               carma_obj<cuFloatComplex> *d_fttotim, long nxsub, long nvalid,
               long npix, long nphase, long nrebin, long nfft, long ntot,
               long npup, float pdiam, float nphotons, float nphot4imat,
               int lgs, bool is_low_order, bool roket, int device);
  sutra_wfs_sh(const sutra_wfs_sh &wfs);
  ~sutra_wfs_sh();

  int define_mpi_rank(int rank, int size);
  int allocate_buffers(map<vector<int>, cufftHandle *> campli_plans,
                       map<vector<int>, cufftHandle *> fttotim_plans);

  int loadarrays(int *phasemap, int *hrmap, int *binmap, float *offsets,
                 float *fluxPerSub, int *validsubsx, int *validsubsy,
                 int *istart, int *jstart, cuFloatComplex *kernel);

  int fill_binimage(int async);
  int comp_image(bool noise = true);

 private:
  int comp_generic();
};

#endif  // _SUTRA_WFS_SH_H_
