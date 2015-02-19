#ifndef _SUTRA_DM_H_
#define _SUTRA_DM_H_

#include <map>
#include <sutra_phase.h>
#include <sutra_kl.h>
#include <sutra_ao_utils.h>
#include <carma_utils.h>

typedef std::pair<std::string, float> type_screen;

using namespace std;

class sutra_dm {
public:
  int device;
  string type;
  long ninflu;
  long influsize;
  long dim;
  float push4imat;

  sutra_phase *d_shape;

  carma_obj<float> *d_comm;

  carma_obj<float> *d_influ; // if relevant

  // pzt 
  carma_obj<int> *d_influpos;
  carma_obj<int> *d_npoints;
  carma_obj<int> *d_istart;
  carma_obj<int> *d_xoff;
  carma_obj<int> *d_yoff;
  carma_obj<float> *d_KLbasis;
  //carma_sparse_obj<float> *d_IFsparse;
  //carma_obj<float> *d_commdouble;
  //carma_obj<double> *d_shapedouble;

  //zernike
  carma_obj<float> *d_coeffs;
  carma_obj<float> *d_mask;
  carma_obj<float> *d_zr;
  carma_obj<float> *d_ztheta;

  sutra_kl *d_kl;

  carma_context *current_context;
  cublasHandle_t cublas_handle() {
      return current_context->get_cublasHandle();
    }
  cusparseHandle_t cusparse_handle() {
        return current_context->get_cusparseHandle();
      }

public:
  sutra_dm(carma_context *context, const char* type, long dim, long ninflu,
      long influsize, long ninflupos, long n_npoints, float push4imat,
      int device);
  sutra_dm(const sutra_dm& dm);
  ~sutra_dm();

  int
  pzt_loadarrays(float *influ, int *influpos, int *npoints, int *istart,
      int *xoff, int *yoff);
  int
  kl_loadarrays(float *rabas, float *azbas, int *ord, float *cr, float *cp);
  int
  reset_shape();
  int
  comp_shape();
  int
  comp_shape(float *comm);
  int
  comp_oneactu(int nactu, float ampli);
  // Florian features
  int
  kl_floloadarrays(float *covmat, float *filter, float *evals, float *bas);

  template<class T>
  int
  get_IF(T *IF, int *indx_pup, long nb_pts, float ampli);
  template<class T>
  int
  get_IF_sparse(carma_sparse_obj<T> *&d_IFsparse, int *indx_pup, long nb_pts, float ampli, int puponly);

  int
  do_geomat(float *d_geocov, float *d_IF, long n_pts);

  template<class T>
  int
  do_geomatFromSparse(T *d_geocov, carma_sparse_obj<T> *d_IFsparse);

  int
  DDiago(carma_obj<float> *d_statcov, carma_obj<float>*d_geocov);
  int
  compute_KLbasis(float *xpos, float *ypos, int *indx, long dim, float norm, float ampli);
  int
  piston_filt(carma_obj <float> *d_statcov);
  int
  set_comkl(float *comvec);
};

class sutra_dms {
public:
  int ndm;
  map<type_screen, sutra_dm *> d_dms;

public:
  sutra_dms(int ndm);
  ~sutra_dms();

  int
  add_dm(carma_context *context, const char* type, float alt, long dim,
      long ninflu, long influsize, long ninflupos, long n_npoints,
      float push4imat, int device);
  int
  remove_dm(const char* type, float alt);

  int nact_total();
};

template<class T>
void
comp_dmshape(int threads, int blocks, T *d_idata, T *d_odata, int *pos,
    int *istart, int *npts, T *comm, unsigned int n, int N);
template<class T>
void
oneactu(int threads, int blocks, T *d_idata, T *d_odata, int nactu, T ampli,
    int *xoff, int *yoff, int dim_im, int dim_influ, int N);
template<class T>
void
oneactu(int threads, int blocks, T *d_idata, T *d_odata, int nactu, T ampli,
    int dim_im, int dim_influ, int N);
template<class T>
void
comp_fulldmshape(int threads, int blocks, T *d_idata, T *d_odata, int ninflu,
    int diminflu, T *comm, int N);

template<class T>
int
getIF(T *IF, float *dmshape, int *indx_pup, long nb_pts, int column, long nb_col, int puponly, carma_device *device);
int
dm_dostatmat(float *d_statcov, long Nkl, float *d_xpos, float *d_ypos, float norm, carma_device *device);
int
fill_filtermat(float *filter, int nactu, int N, carma_device *device);
int
find_nnz(float *d_data, int N,carma_device *device);

#endif // _SUTRA_DM_H_
