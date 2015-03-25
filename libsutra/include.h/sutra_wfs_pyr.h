#ifndef _SUTRA_WFS_PYR_H_
#define _SUTRA_WFS_PYR_H_

#include <vector>
#include <map>
#include <sutra_telemetry.h>
#include <sutra_target.h>
#include <sutra_phase.h>
#include <sutra_lgs.h>
#include <sutra_wfs.h>

using namespace std;

class sutra_wfs_pyr : public sutra_wfs {
public:
  // pyramid only
  carma_obj<float> *d_hrimg;
  carma_obj<float> *d_submask;
  carma_obj<float> *d_psum;
  carma_obj<cuFloatComplex> *d_phalfxy;
  carma_obj<cuFloatComplex> *d_poffsets;

  carma_host_obj<int> *pyr_cx;
  carma_host_obj<int> *pyr_cy;

public:
  sutra_wfs_pyr(carma_context *context,sutra_sensors *sensors, long nxsub, long nvalid,
      long npix, long nphase, long nrebin, long nfft, long ntot, long npup,
      float pdiam, float nphotons, int lgs, int device);
  sutra_wfs_pyr(const sutra_wfs_pyr& wfs);
  ~sutra_wfs_pyr();

  int
  wfs_initarrays(cuFloatComplex *halfxy, cuFloatComplex *offsets,
      float *focmask, float *pupil, int *isvalid, int *cx, int *cy,
      float *sincar, int *phasemap, int *validsubsx, int *validsubsy);

  virtual int
  comp_image()=0;

  int define_mpi_rank(int rank, int size) {return EXIT_SUCCESS;}
  int allocate_buffers(sutra_sensors *sensors) {return EXIT_SUCCESS;}

private:
  virtual int
  comp_generic()=0;
};

#endif // _SUTRA_WFS_PYR_H_
