#ifndef SUTRA_CONTROLLER_GEO_H_
#define SUTRA_CONTROLLER_GEO_H_

#include <sutra_controller.h>

template <typename T>
class sutra_controller_geo : public sutra_controller<T> {
 public:
  T gain;
  long Nphi;
  int Ntt;

  carma_obj<T> *d_gain;
  carma_obj<T> *d_proj;
  carma_obj<double> *d_phi;
  carma_obj<T> *d_phif;
  carma_obj<int> *d_indx_pup;
  carma_obj<int> *d_indx_mpup;
  carma_sparse_obj<double> *d_IFsparse;
  carma_obj<T> *d_geocov;
  carma_obj<double> *d_compdouble;
  carma_obj<float> *d_compfloat;
  carma_obj<T> *d_TT;
  carma_obj<T> *d_geocovTT;
  //  carma_obj<T> *d_Btt;
  // carma_obj<T> *d_cenbuff; // centroids circular buffer

 public:
  sutra_controller_geo(carma_context *context, long nactu, long Nphi, T delay,
                       sutra_dms *dms, int *idx_dms, int ndm,
                       bool wfs_direction);
  sutra_controller_geo(const sutra_controller_geo &controller);
  ~sutra_controller_geo();

  string get_type();

  cusparseHandle_t cusparse_handle() {
    return this->current_context->get_cusparseHandle();
  }
  int load_Btt(T *Btt_pzt, T *Btt_TT);
  int set_gain(T gain);
  int load_mgain(T *mgain);
  int comp_dphi(sutra_source *target, bool wfs_direction);
  int comp_com();
  int init_proj(sutra_dms *dms, int *indx_dm, T *unitpervolt, int *indx_pup);
  int init_proj_sparse(sutra_dms *dms, int *indx_dm, T *unitpervolt,
                       int *indx_pup, int *indx_mpup, bool roket);
};

#endif /* SUTRA_CONTROLLER_GEO_H_ */
