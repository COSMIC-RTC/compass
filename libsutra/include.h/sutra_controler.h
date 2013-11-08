#ifndef _SUTRA_CONTROLER_H_
#define _SUTRA_CONTROLER_H_

#include <carma_cublas.h>
#include <carma_host_obj.h>
#include <sutra_wfs.h>
#include <sutra_dm.h>
#include <sutra_centroider.h>
#include <sutra_ao_utils.h>

using namespace std;

class sutra_controler {
 public:
  string                   typec;
  int                      device;
  int                      nvalid;
  int                      nactu;
  int                      delay;
  float                    gain;

  carma_obj<float>          *d_imat;     
  carma_obj<float>          *d_cmat;     
  carma_obj<float>          *d_gain;     

  // svd computations
  carma_host_obj<float>     *h_imat;     
  carma_host_obj<float>     *h_cmat;
  carma_obj<float>     	   *d_eigenvals;
  carma_host_obj<float>     *h_eigenvals;
  carma_host_obj<float>     *h_mateigens;
  carma_obj<float>          *d_mateigens;
  carma_host_obj<float>     *h_mod2act;
  carma_obj<float>          *d_mod2act;
  carma_host_obj<float>     *h_mes2mod;
  carma_obj<float>          *d_mes2mod;
  carma_obj<float>          *d_tmp;

  // loop components
  carma_obj<float>          *d_centroids; // current centroids    
  carma_obj<float>          *d_cenbuff;   // centroids circular buffer  
  carma_obj<float>          *d_com;       // current command
  //carma_obj<float>          *d_combuff;   // command circular buffer 
  carma_obj<float>          *d_err;       // current error
  //carma_obj<float>          *d_err;       // error circular buffer  

  carma_context *current_context;

  carma_streams 		    *streams;
  int  			    nstreams;
  cublasHandle_t            cublas_handle;

 public:
  sutra_controler(carma_context *context, long nvalid, long nactu, long delay, int device, const char* typec);
  sutra_controler(const sutra_controler& controler);
  ~sutra_controler();

  int svdec_imat();
  int build_cmat(int nfilt,bool filt_tt);
  int build_cmat(int nfilt);
  int frame_delay();
  int comp_com();
  int set_gain(float gain);
  int load_mgain(float *mgain);
  int set_delay(int delay);

};

int shift_buf(float *d_data,int offset, int N,int device);
int mult_vect(float *d_data,float *scale, int N,int device);
int mult_int(float *o_data,float *i_data,float *scale, float gain, int N,int device);
int mult_int(float *o_data,float *i_data,float *scale, float gain, int N,int device, carma_streams *streams);

#endif // _SUTRA_CONTROLER_H_

