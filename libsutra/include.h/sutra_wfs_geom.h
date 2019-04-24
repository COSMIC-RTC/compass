/**
 * \file sutra_wfs_atmos.h
 *
 * \class sutra_wfs_atmos
 *
 * \ingroup libsutra
 *
 * \brief this class provides the wfs_atmos features to COMPASS
 *
 * \authors Damien Gratadour & Arnaud Sevin & Florian Ferreira
 *
 * \version 1.0
 *
 * \date 2011/01/28
 *
 */
#ifndef _SUTRA_WFS_GEOM_H_
#define _SUTRA_WFS_GEOM_H_

#include <sutra_lgs.h>
#include <sutra_phase.h>
#include <sutra_target.h>
#include <sutra_telemetry.h>
#include <sutra_telescope.h>
#include <sutra_wfs.h>
#include <map>
#include <vector>

using std::string;
class sutra_wfs_geom : public sutra_wfs {
 public:
 public:
  sutra_wfs_geom(carma_context *context, sutra_telescope *d_tel, long nxsub,
                 long nvalid, long nphase, long npup, float pdiam, int device);
  sutra_wfs_geom(const sutra_wfs_geom &wfs);
  ~sutra_wfs_geom();

  int wfs_initarrays(int *phasemap, float *offsets, float *fluxPerSub,
                     int *validsubsx, int *validsubsy);
  int slopes_geom(int type, float *slopes);
  int slopes_geom(int type);

  int define_mpi_rank(int rank, int size) { return EXIT_SUCCESS; }
  int allocate_buffers(map<vector<int>, cufftHandle *> campli_plans,
                       map<vector<int>, cufftHandle *> fttotim_plans) {
    return EXIT_SUCCESS;
  }

  int fill_binimage(int async) { return 0; }
  int comp_image() { return 0; }

 private:
  int comp_generic() { return 0; }
};

#endif  // _SUTRA_WFS_GEOM_H_
