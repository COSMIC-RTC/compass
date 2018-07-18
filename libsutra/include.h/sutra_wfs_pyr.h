#ifndef _SUTRA_WFS_PYR_H_
#define _SUTRA_WFS_PYR_H_

#include <sutra_lgs.h>
#include <sutra_phase.h>
#include <sutra_target.h>
#include <sutra_telemetry.h>
#include <sutra_telescope.h>
#include <sutra_wfs.h>
#include <map>
#include <vector>

using std::string;

class sutra_wfs_pyr : public sutra_wfs {
 public:
  // pyramid only
  carma_obj<float> *d_hrimg;
  carma_obj<float> *d_submask;
  carma_obj<float> *d_psum;
  carma_obj<cuFloatComplex> *d_phalfxy;
  carma_obj<cuFloatComplex> *d_poffsets;

  carma_host_obj<float> *pyr_cx;
  carma_host_obj<float> *pyr_cy;

 public:
  sutra_wfs_pyr(carma_context *context, sutra_telescope *d_tel,
                carma_obj<cuFloatComplex> *d_camplipup,
                carma_obj<cuFloatComplex> *d_camplifoc,
                carma_obj<cuFloatComplex> *d_fttotim, long nxsub, long nvalid,
                long npix, long nphase, long nrebin, long nfft, long ntot,
                long npup, float pdiam, float nphotons, float nphot4imat,
                int lgs, bool roket, int device, const char *type_pyr);
  sutra_wfs_pyr(const sutra_wfs_pyr &wfs);
  ~sutra_wfs_pyr();

  int loadarrays(cuFloatComplex *halfxy, cuFloatComplex *offsets,
                 float *focmask, float *cx, float *cy, float *sincar,
                 int *phasemap, int *validsubsx, int *validsubsy);

  int fill_binimage(int async);
  virtual int comp_image() = 0;

  int define_mpi_rank(int rank, int size) { return EXIT_SUCCESS; }
  int allocate_buffers(map<vector<int>, cufftHandle *> campli_plans,
                       map<vector<int>, cufftHandle *> fttotim_plans) {
    return EXIT_SUCCESS;
  }

 private:
  virtual int comp_generic() = 0;
};

#endif  // _SUTRA_WFS_PYR_H_
