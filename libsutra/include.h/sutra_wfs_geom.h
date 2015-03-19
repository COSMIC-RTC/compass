#ifndef _SUTRA_WFS_GEOM_H_
#define _SUTRA_WFS_GEOM_H_

#include <vector>
#include <map>
#include <sutra_telemetry.h>
#include <sutra_target.h>
#include <sutra_phase.h>
#include <sutra_lgs.h>
#include <sutra_wfs.h>

using namespace std;
class sutra_sensors;
class sutra_wfs_geom : public sutra_wfs{
public:


public:
    sutra_wfs_geom(carma_context *context, long nxsub, long nvalid, long nphase,
      long npup, float pdiam, int device);
    sutra_wfs_geom(const sutra_wfs_geom& wfs);
  ~sutra_wfs_geom();

  int
  wfs_initarrays(int *phasemap, float *offsets, float *pupil, float *fluxPerSub,
      int *isvalid, int *validsubsx, int *validsubsy);
  int
  slopes_geom(int type, float *slopes);
  int
  slopes_geom(int type);

  int
  comp_image(){return 0;}
  int
  comp_image_tele(){return 0;}

private:
  int
  comp_generic(){return 0;}
};


#endif // _SUTRA_WFS_GEOM_H_