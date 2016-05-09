#ifndef _SUTRA_WFS_PYR_PYR4_H_
#define _SUTRA_WFS_PYR_PYR4_H_

#include <vector>
#include <map>
#include <sutra_telemetry.h>
#include <sutra_target.h>
#include <sutra_phase.h>
#include <sutra_lgs.h>
#include <sutra_wfs_pyr.h>
#include <sutra_telescope.h>

using namespace std;
class sutra_sensors;
class sutra_wfs_pyr_pyr4: public sutra_wfs_pyr {
  public:

  public:
    sutra_wfs_pyr_pyr4(carma_context *context, sutra_telescope *d_tel,
                       sutra_sensors *sensors, long nxsub, long nvalid,
                       long npix, long nphase, long nrebin, long nfft,
                       long ntot, long npup, float pdiam, float nphotons,
                       float nphot4imat, int lgs, int device);
    sutra_wfs_pyr_pyr4(const sutra_wfs_pyr_pyr4& wfs);
    ~sutra_wfs_pyr_pyr4();

    int comp_image();

  private:
    int comp_generic();
};

#endif // _SUTRA_WFS_PYR_PYR4_H_
