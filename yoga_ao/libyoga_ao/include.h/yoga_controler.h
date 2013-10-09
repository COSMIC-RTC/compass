#ifndef _YOGA_CONTROLER_H_
#define _YOGA_CONTROLER_H_

#include <yoga_cublas.h>
#include <yoga_host_obj.h>
#include <yoga_wfs.h>
#include <yoga_dm.h>
#include <yoga_centroider.h>
#include <yoga_ao_utils.h>

using namespace std;

class yoga_controler {
 public:
  string                   typec;
  int                      device;
  int                      nvalid;
  int                      nactu;
  int                      delay;
  float                    gain;

  yoga_obj<float>          *d_imat;     
  yoga_obj<float>          *d_cmat;     
  yoga_obj<float>          *d_gain;     

  // svd computations
  yoga_host_obj<float>     *h_imat;     
  yoga_host_obj<float>     *h_cmat;
  yoga_obj<float>     	   *d_eigenvals;
  yoga_host_obj<float>     *h_eigenvals;
  yoga_host_obj<float>     *h_mateigens;
  yoga_obj<float>          *d_mateigens;
  yoga_host_obj<float>     *h_mod2act;
  yoga_obj<float>          *d_mod2act;
  yoga_host_obj<float>     *h_mes2mod;
  yoga_obj<float>          *d_mes2mod;
  yoga_obj<float>          *d_tmp;

  // loop components
  yoga_obj<float>          *d_centroids; // current centroids    
  yoga_obj<float>          *d_cenbuff;   // centroids circular buffer  
  yoga_obj<float>          *d_com;       // current command
  //yoga_obj<float>          *d_combuff;   // command circular buffer 
  yoga_obj<float>          *d_err;       // current error
  //yoga_obj<float>          *d_err;       // error circular buffer  

  yoga_context *current_context;

  yoga_streams 		    *streams;
  int  			    nstreams;
  cublasHandle_t            cublas_handle;

 public:
  yoga_controler(yoga_context *context, long nvalid, long nactu, long delay, int device, const char* typec);
  yoga_controler(const yoga_controler& controler);
  ~yoga_controler();

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
int mult_int(float *o_data,float *i_data,float *scale, float gain, int N,int device, yoga_streams *streams);

#endif // _YOGA_CONTROLER_H_

