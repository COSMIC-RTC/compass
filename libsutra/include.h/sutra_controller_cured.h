/**
 * \file sutra_controller_cured.h
 *
 * \class sutra_controller_cured
 *
 * \ingroup libsutra
 *
 * \brief this class provides the controller_cured features to COMPASS
 *
 * \authors Damien Gratadour & Arnaud Sevin & Florian Ferreira
 *
 * \version 4.2.0
 *
 * \date 2011/01/28
 *
 */
#ifndef _sutra_controller_cured_H_
#define _sutra_controller_cured_H_

#include <sutra_controller.h>

template <typename Tcomp, typename Tout>
class sutra_controller_cured : public sutra_controller<Tcomp, Tout> {
 public:
  int ndivs;     // number of subdivision levels for cured
  bool tt_flag;  // flag for separate tt

  // data for CuReD */
  carma_host_obj<Tcomp> *h_centroids;
  carma_host_obj<Tcomp> *h_err;
  carma_obj<Tcomp> *d_err;      // current error
  carma_obj<Tcomp> *d_cenbuff;  // centroids circular buffer

  // data for CuReD */
  carma_obj<Tcomp> *d_imat;

  // structures needed to run CuReD */
  // sysCure* h_syscure;
  void *h_syscure;
  // parCure* h_parcure;
  void *h_parcure;

 public:
  sutra_controller_cured(carma_context *context, long nvalid, long nslope,
                         long nactu, float delay, sutra_dms *dms, int *idx_dms,
                         int ndm);
  sutra_controller_cured(const sutra_controller_cured &controller);
  ~sutra_controller_cured();

  string get_type() { return "cured"; }

  int comp_com();

  int init_cured(int nxsubs, int *isvalid, int ndivs, int tt);
  int frame_delay();
};

#endif  // _sutra_controller_cured_H_
