#ifndef _SUTRA_GROOT_H_
#define _SUTRA_GROOT_H_

#include <carma.h>
#include <carma_obj.h>
#include <carma_sparse_obj.h>
#include <carma_host_obj.h>
#include <carma_cublas.h>

template<class T_data>
class sutra_groot {
  public:
    carma_context *current_context;
    int device;

    int nactus; // number of actuators
    int nlayers; // number of atmos layers
    T_data gsangle; // Guide star angle [rad]
    T_data fc; // DM cut-off frequency [m]

    carma_host_obj<T_data> *h_vdt; // v*dt/g [m/s]
    carma_host_obj<T_data> *h_Htheta; // H*theta (theta = GS radial distance) [m]
    carma_host_obj<T_data> *h_L0; // L0 [m]
    carma_host_obj<T_data> *h_winddir; // wind directions [rad]
    carma_host_obj<T_data> *h_scale; // r0**(-5/3) * frac * (lambda/2pi)**2

    //Covariance matrix estimation
    carma_obj<T_data> *d_Cerr; // Residual error covariance matrix on DM actuators
    carma_obj<T_data> *d_TT; // TT component of the residual error covariance matrix

    carma_obj<T_data> *d_TTPfilter; // Tip-tilt and piston filter matrix (= Btt.dot(P))
    carma_obj<T_data> *d_pzt2tt; // pzt to TT matrix
    carma_obj<T_data> *d_Nact; // Coupling matrix
    carma_obj<T_data> *d_xpos; // X-positions of DM actuators [m]
    carma_obj<T_data> *d_ypos; // Y-positions of DM actuators [m]

    // Dphi lowpass
    carma_obj<T_data> *d_tab_int_x; // Tabulated integral
    carma_obj<T_data> *d_tab_int_y;


  public:
    sutra_groot(carma_context *context, int device, int nactus,
                int nlayers, T_data gsangle, T_data *vdt,
                T_data *Htheta, T_data *L0, T_data *winddir, T_data *scale,
                T_data *pzt2tt, T_data *TTPfilter,T_data *Nact,
                T_data *xpos, T_data *ypos, T_data fc);
    ~sutra_groot();

    cublasHandle_t cublas_handle() {
      return current_context->get_cublasHandle();
    }

    int compute_Cerr();
};
template<class T_data>
int compute_Cerr_layer(T_data *Cerr, int N, T_data *tab_int_x, T_data *tab_int_y,
                            T_data *xpos, T_data *ypos, T_data vdt,
                            T_data Htheta, T_data L0, T_data fc,
                            T_data winddir, T_data gsangle,
                            T_data scale, int Ntab, carma_device *device);
template<class T_data>
int tab_u831J0(T_data *tab_int_x, T_data *tab_int_y, int N, carma_device *device);
template<class T_data>
void cumsum(T_data *odata, T_data *idata, int N);
template<class T_data>
int add_transpose(T_data *Cerr, int N, carma_device *device);

#endif
