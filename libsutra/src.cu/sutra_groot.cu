// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.
//  Distributed under GNU - LGPL
//
//  COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
//  General Public License as published by the Free Software Foundation, either version 3 of the License,
//  or any later version.
//
//  COMPASS: End-to-end AO simulation tool using GPU acceleration
//  The COMPASS platform was designed to meet the need of high-performance for the simulation of AO systems.
//
//  The final product includes a software package for simulating all the critical subcomponents of AO,
//  particularly in the context of the ELT and a real-time core based on several control approaches,
//  with performances consistent with its integration into an instrument. Taking advantage of the specific
//  hardware architecture of the GPU, the COMPASS tool allows to achieve adequate execution speeds to
//  conduct large simulation campaigns called to the ELT.
//
//  The COMPASS platform can be used to carry a wide variety of simulations to both testspecific components
//  of AO of the E-ELT (such as wavefront analysis device with a pyramid or elongated Laser star), and
//  various systems configurations such as multi-conjugate AO.
//
//  COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
//  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//  See the GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License along with COMPASS.
//  If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>.
// -----------------------------------------------------------------------------

//! \file      sutra_groot.cu
//! \ingroup   libsutra
//! \class     SutraGroot
//! \brief     this class provides the groot features to COMPASS
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.3.0
//! \date      2022/01/24
//! \copyright GNU Lesser General Public License

#include <sutra_groot.h>

template <class T_data, typename Fn>
__device__ T_data macdo_x56_gpu_gb(Fn const &ptr_pow, T_data x, int k)
/* DOCUMENT  macdo_x56_gpu_gb(x)

 Computation of the function
 f(x) = x^(5/6)*K_{5/6}(x)
 using a series for the esimation of K_{5/6}, taken from Rod Conan thesis :
 K_a(x)=1/2 \sum_{n=0}^\infty \frac{(-1)^n}{n!}
 \left(\Gamma(-n-a) (x/2)^{2n+a} + \Gamma(-n+a) (x/2)^{2n-a} \right) ,
 with a = 5/6.

 Setting x22 = (x/2)^2, setting uda = (1/2)^a, and multiplying by x^a,
 this becomes :
 x^a * Ka(x) = 0.5 $ -1^n / n! [ G(-n-a).uda x22^(n+a) + G(-n+a)/uda x22^n ]
 Then we use the following recurrence formulae on the following quantities :
 G(-(n+1)-a) = G(-n-a) / -a-n-1
 G(-(n+1)+a) = G(-n+a) /  a-n-1
 (n+1)! = n! * (n+1)
 x22^(n+1) = x22^n * x22
 and at each iteration on n, one will use the values already computed at step
 (n-1). The values of G(a) and G(-a) are hardcoded instead of being computed.

 The first term of the series has also been skipped, as it
 vanishes with another term in the expression of Dphi.

 SEE ALSO:
 */
{
  const T_data a = 5. / 6.;
  const T_data x2a = ptr_pow(x, (T_data)2. * a), x22 = x * x / 4.;
  T_data x2n;  // x^2.a, etc
  T_data s = 0.0;
  int n;

  const T_data Ga[11] = {0,
                         12.067619015983075,
                         5.17183672113560444,
                         0.795667187867016068,
                         0.0628158306210802181,
                         0.00301515986981185091,
                         9.72632216068338833e-05,
                         2.25320204494595251e-06,
                         3.93000356676612095e-08,
                         5.34694362825451923e-10,
                         5.83302941264329804e-12};

  const T_data Gma[11] = {-3.74878707653729304,     -2.04479295083852408,
                          -0.360845814853857083,    -0.0313778969438136685,
                          -0.001622994669507603,    -5.56455315259749673e-05,
                          -1.35720808599938951e-06, -2.47515152461894642e-08,
                          -3.50257291219662472e-10, -3.95770950530691961e-12,
                          -3.65327031259100284e-14};

  x2n = 0.5;  // init (1/2) * x^0

  s = Gma[0] * x2a;
  s *= x2n;

  // prepare recurrence iteration for next step
  x2n *= x22;  // x^n

#pragma unroll
  for (n = 1; n <= 10; n++) {
    s += (Gma[n] * x2a + Ga[n]) * x2n;
    // prepare recurrence iteration for next step
    x2n *= x22;  // x^n
  }
  return s;
}

__device__ float macdo_x56_gpu(float x, int k) {
  return macdo_x56_gpu_gb<float>(powf, x, k);
}

__device__ double macdo_x56_gpu(double x, int k) {
  return macdo_x56_gpu_gb<double, double (*)(double, double)>(pow, x, k);
}

//------------------------------------------------------------------------------------
template <class T_data, typename Fn, typename Fne>
__device__ T_data asymp_macdo_gpu_gb(Fn const &ptr_pow, Fne const &ptr_exp,
                                     T_data x)
/* DOCUMENT asymp_macdo_gpu_gb(x)

 Computes a term involved in the computation of the phase struct
 function with a finite outer scale according to the Von-Karman
 model. The term involves the MacDonald function (modified bessel
 function of second kind) K_{5/6}(x), and the algorithm uses the
 asymptotic form for x ~ infinity.
 Warnings :
 - This function makes a T_dataing point interrupt for x=0
 and should not be used in this case.
 - Works only for x>0.

 SEE ALSO:
 */
{
  // k2 is the value for
  // gamma_R(5./6)*2^(-1./6)
  const T_data k2 = 1.00563491799858928388289314170833;
  const T_data k3 = 1.25331413731550012081;   //  sqrt(pi/2)
  const T_data a1 = 0.22222222222222222222;   //  2/9
  const T_data a2 = -0.08641975308641974829;  //  -7/89
  const T_data a3 = 0.08001828989483310284;   // 175/2187
  T_data res;
  T_data x_1;

  x_1 = 1. / x;
  res = k2 - k3 * ptr_exp(-x) * ptr_pow(x, (T_data)(1 / 3.)) *
                 (1.0 + x_1 * (a1 + x_1 * (a2 + x_1 * a3)));
  return res;
}

__device__ float asymp_macdo_gpu(float x) {
  return asymp_macdo_gpu_gb<float>(powf, expf, x);
}

__device__ double asymp_macdo_gpu(double x) {
  return asymp_macdo_gpu_gb<double, double (*)(double, double),
                            double (*)(double)>(pow, exp, x);
}

//------------------------------------------------------------------------------------
template <class T_data, typename Fn>
__device__ T_data rodconan_gpu_gb(Fn const &ptr_pow, T_data r, T_data L0, int k)
/* DOCUMENT rodconan_gpu_gb(r,L0,k=)
 The phase structure function is computed from the expression
 Dphi(r) = k1  * L0^(5./3) * (k2 - (2.pi.r/L0)^5/6 K_{5/6}(2.pi.r/L0))

 For small r, the expression is computed from a development of
 K_5/6 near 0. The value of k2 is not used, as this same value
 appears in the series and cancels with k2.
 For large r, the expression is taken from an asymptotic form.

 SEE ALSO:
 */
{
  T_data res = 0;

  // k1 is the value of :
  // 2*gamma_R(11./6)*2^(-5./6)*pi^(-8./3)*(24*gamma_R(6./5)/5.)^(5./6);
  const T_data k1 = 0.1716613621245709486;
  const T_data dprf0 = (2 * CARMA_PI / L0) * r;
  // k2 is the value for gamma_R(5./6)*2^(-1./6),
  // but is now unused
  // k2 = 1.0056349179985892838;

  // Xlim = 0.75*2*pi;   // = 4.71239
  if (dprf0 > 4.71239)
    res = asymp_macdo_gpu(dprf0);
  else
    res = -macdo_x56_gpu(dprf0, k);

  res *= k1 * ptr_pow(L0, (T_data)5. / 3);

  return res;
}

__device__ float rodconan_gpu(float r, float L0, int k) {
  return rodconan_gpu_gb<float>(powf, r, L0, k);
}

__device__ double rodconan_gpu(double r, double L0, int k) {
  return rodconan_gpu_gb<double, double (*)(double, double)>(pow, r, L0, k);
}

template <class T_data>
__device__ T_data unMoinsJ0(T_data x) {
  if (x < 0.1) {
    T_data x22 = (x * x) / 4.;
    return (1.0 - x22 / 4.) * x22;
  } else
    return 1.0f - j0(x);
}

template __device__ float unMoinsJ0<float>(float x);
template __device__ double unMoinsJ0<double>(double x);

template <class T_data, typename Fn>
__device__ void compute_u831J0_gen(Fn const &ptr_exp, T_data *x, T_data *y,
                                   int npts, T_data tmin, T_data tmax,
                                   T_data dt) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  T_data t;
  while (tid < npts) {
    t = tmin + tid * dt;
    x[tid] = ptr_exp(t);
    y[tid] = ptr_exp(-t * (5 / 3.)) * unMoinsJ0(ptr_exp(t)) * dt;
    tid += blockDim.x * gridDim.x;
  }
}

template <class T_data>
__global__ void compute_u831J0(T_data *x, T_data *y, int npts, T_data tmin,
                               T_data tmax, T_data dt, CarmaDevice *device);
template <>
__global__ void compute_u831J0<float>(float *x, float *y, int npts, float tmin,
                                      float tmax, float dt,
                                      CarmaDevice *device) {
  compute_u831J0_gen<float>(expf, x, y, npts, tmin, tmax, dt);
}
template <>
__global__ void compute_u831J0<double>(double *x, double *y, int npts,
                                       double tmin, double tmax, double dt,
                                       CarmaDevice *device) {
  compute_u831J0_gen<double, double (*)(double)>(exp, x, y, npts, tmin, tmax,
                                                 dt);
}

template <class T_data>
__global__ void intfrominftomin(T_data *data, T_data smallInt, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < N) {
    data[tid] += smallInt;
    tid += blockDim.x * gridDim.x;
  }
}

template __global__ void intfrominftomin<float>(float *data, float smallInt,
                                                int N);
template __global__ void intfrominftomin<double>(double *data, double smallInt,
                                                 int N);

template <class T_data>
__global__ void cuda_zcen_krnl(T_data *idata, T_data *odata, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    odata[tid] = (idata[tid + 1] + idata[tid]) / 2.0;
    tid += blockDim.x * gridDim.x;
  }
}

template __global__ void cuda_zcen_krnl<float>(float *idata, float *odata,
                                               int N);
template __global__ void cuda_zcen_krnl<double>(double *idata, double *odata,
                                                int N);

template <class T_data>
void cumsum(T_data *odata, T_data *idata, int N) {
  odata[0] = 0;
  for (int i = 1; i < N; i++) {
    odata[i] = idata[i - 1] + odata[i - 1];
  }
}

template void cumsum<float>(float *odata, float *idata, int N);
template void cumsum<double>(double *odata, double *idata, int N);

template <class T_data>
int tab_u831J0(T_data *tab_int_x, T_data *tab_int_y, int npts,
               CarmaDevice *device) {
  T_data tmin = -4.;
  T_data tmax = 10.;
  T_data *t = (T_data *)malloc(sizeof(T_data) * npts);
  cudaMemcpy(t, tab_int_x, (npts) * sizeof(T_data), cudaMemcpyDeviceToHost);
  carma_check_msg("copy test");

  T_data *temp;
  temp = (T_data *)malloc((npts - 1) * sizeof(T_data));
  T_data *tab;
  tab = (T_data *)malloc((npts) * sizeof(T_data));
  T_data *temp_d;
  T_data dt = (tmax - tmin) / (npts - 1);

  int nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, npts, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  compute_u831J0<<<grid, threads>>>(tab_int_x, tab_int_y, npts, tmin, tmax, dt,
                                    device);
  carma_check_msg("compute_u831J0<<<>>> execution failed\n");
  // DEBUG_TRACE("tab_int !\n");
  cudaMalloc((void **)&(temp_d), (npts - 1) * sizeof(T_data));
  carma_check_msg("alloc gpu temp_d");

  get_num_blocks_and_threads(device, npts - 1, nb_blocks, nb_threads);
  dim3 grid2(nb_blocks), threads2(nb_threads);
  cuda_zcen_krnl<<<grid2, threads2>>>(tab_int_y, temp_d, npts - 1);
  carma_check_msg("cuda_zcen_krnl<<<>>> execution failed\n");
  // cuda_zcen(tab_int_y,temp_d, npts-1, device);
  // DEBUG_TRACE("tab_int !\n");
  cudaMemcpy(temp, temp_d, (npts - 1) * sizeof(T_data), cudaMemcpyDeviceToHost);
  carma_check_msg("copy cpu temp");
  cudaFree(temp_d);
  cumsum(tab, temp, npts);
  // DEBUG_TRACE("tab_int !\n");
  cudaMemcpy(tab_int_y, tab, (npts) * sizeof(T_data), cudaMemcpyHostToDevice);
  carma_check_msg("copy gpu tab");
  // DEBUG_TRACE("tab_int !\n");
  T_data smallx = exp(tmin);
  T_data smallInt =
      (0.75 * pow(smallx, (T_data)(1 / 3.)) * (1 - smallx * smallx / 112.));
  // DEBUG_TRACE("tab_int !\n");

  intfrominftomin<<<grid, threads>>>(tab_int_y, smallInt, npts);
  carma_check_msg("intfrominftomin<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}
template int tab_u831J0(float *tab_int_x, float *tab_int_y, int npts,
                        CarmaDevice *device);

template int tab_u831J0(double *tab_int_x, double *tab_int_y, int npts,
                        CarmaDevice *device);

template <class T_data, typename Fn, typename Fne, typename Fnl>
__device__ T_data Ij0t83_gen(Fn const &ptr_pow, Fne const &ptr_exp,
                             Fnl const &ptr_log, T_data x, T_data *tab_x,
                             T_data *tab_y, long npts) {
  if (x <= ptr_exp(-3.0))
    return (T_data)(0.75 * ptr_pow((T_data)x, 1 / 3.) * (1 - x * x / 112.));
  else {
    T_data dt = 14. / (npts - 1);  // 14 = tmax - tmin
    T_data convert = (ptr_log(x) + 4.) / dt;
    long i0 = (long)convert;
    long i1 = i0 + 1;

    return (((x - tab_x[i0]) * tab_y[i1] + (tab_x[i1] - x) * tab_y[i0])) /
           (tab_x[i1] - tab_x[i0]);
  }
}

__device__ float Ij0t83(float x, float *tab_x, float *tab_y, long npts) {
  return Ij0t83_gen<float>(powf, expf, logf, x, tab_x, tab_y, npts);
}
__device__ double Ij0t83(double x, double *tab_x, double *tab_y, long npts) {
  return Ij0t83_gen<double, double (*)(double, double), double (*)(double),
                    double (*)(double)>(pow, exp, log, x, tab_x, tab_y, npts);
}

template <class T_data, typename Fn>
__device__ T_data DPHI_highpass_gen(Fn const &ptr_pow, T_data r, T_data fc,
                                    T_data *tab_x, T_data *tab_y, long npts) {
  return ptr_pow(r, 5 / 3.) *
         (1.1183343328701949 -
          Ij0t83(2 * CARMA_PI * fc * r, tab_x, tab_y, npts)) *
         ptr_pow(2 * CARMA_PI, 8 / 3.) * 2 * 0.0228956;
}

__device__ float DPHI_highpass(float r, float fc, float *tab_x, float *tab_y,
                               long npts) {
  return DPHI_highpass_gen<float>(powf, r, fc, tab_x, tab_y, npts);
}
__device__ double DPHI_highpass(double r, double fc, double *tab_x,
                                double *tab_y, long npts) {
  return DPHI_highpass_gen<double, double (*)(double, double)>(
      pow, r, fc, tab_x, tab_y, npts);
}

template <class T_data, typename Fn>
__device__ T_data DPHI_highpass_gen(Fn const &ptr_sqrt, T_data x, T_data y,
                                    T_data fc, T_data *tab_x, T_data *tab_y,
                                    long npts) {
  return DPHI_highpass(ptr_sqrt(x * x + y * y), fc, tab_x, tab_y, npts);
}

__device__ float DPHI_highpass(float x, float y, float fc, float *tab_x,
                               float *tab_y, long npts) {
  return DPHI_highpass_gen<float>(sqrtf, x, y, fc, tab_x, tab_y, npts);
}
__device__ double DPHI_highpass(double x, double y, double fc, double *tab_x,
                                double *tab_y, long npts) {
  return DPHI_highpass_gen<double, double (*)(double)>(sqrt, x, y, fc, tab_x,
                                                       tab_y, npts);
}

template <class T_data, typename Fn>
__device__ T_data DPHI_lowpass_gen(Fn const &ptr_sqrt, T_data x, T_data y,
                                   T_data L0, T_data fc, T_data *tab_int_x,
                                   T_data *tab_int_y, long npts) {
  T_data r = ptr_sqrt(x * x + y * y);

  return rodconan_gpu(r, L0, 10) -
         DPHI_highpass(r, fc, tab_int_x, tab_int_y, npts);
}

__device__ float DPHI_lowpass(float x, float y, float L0, float fc,
                              float *tab_int_x, float *tab_int_y, long npts) {
  return DPHI_lowpass_gen<float>(sqrtf, x, y, L0, fc, tab_int_x, tab_int_y,
                                 npts);
}
__device__ double DPHI_lowpass(double x, double y, double L0, double fc,
                               double *tab_int_x, double *tab_int_y,
                               long npts) {
  return DPHI_lowpass_gen<double, double (*)(double)>(
      sqrt, x, y, L0, fc, tab_int_x, tab_int_y, npts);
}

template <class T_data>
__global__ void compute_Ca_element_XX(T_data *CaXX, int N, T_data *tab_int_x,
                                      T_data *tab_int_y, T_data *xpos,
                                      T_data *ypos, T_data d, T_data fc,
                                      T_data scale, T_data weight,
                                      T_data offset, int Ntab) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < N * N) {
    int i = tid / N;
    int j = tid - i * N;
    T_data yoff = offset;
    // if(j < i)
    //   yoff = -offset;
    T_data xij = xpos[i] - xpos[j];
    T_data yij = ypos[i] - ypos[j] + yoff;

    CaXX[tid] +=
        ((DPHI_highpass(xij - d, yij, fc, tab_int_x, tab_int_y, Ntab) +
          DPHI_highpass(xij + d, yij, fc, tab_int_x, tab_int_y, Ntab) -
          2 * DPHI_highpass(xij, yij, fc, tab_int_x, tab_int_y, Ntab)) *
         scale * weight);

    tid += blockDim.x * gridDim.x;
  }
}

template __global__ void compute_Ca_element_XX(float *CaXX, int nssp,
                                               float *tab_int_x,
                                               float *tab_int_y, float *xpos,
                                               float *ypos, float d, float fc,
                                               float scale, float weight,
                                               float offset, int Ntab);
template __global__ void compute_Ca_element_XX(
    double *CaXX, int nssp, double *tab_int_x, double *tab_int_y, double *xpos,
    double *ypos, double d, double fc, double scale, double weight,
    double offset, int Ntab);

template <class T_data>
__global__ void compute_Ca_element_YY(T_data *CaYY, int N, T_data *tab_int_x,
                                      T_data *tab_int_y, T_data *xpos,
                                      T_data *ypos, T_data d, T_data fc,
                                      T_data scale, T_data weight,
                                      T_data offset, int Ntab) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < N * N) {
    int i = tid / N;
    int j = tid - i * N;
    T_data xoff = offset;
    // if(j < i)
    //   xoff = -offset;
    T_data xij = xpos[i] - xpos[j] + xoff;
    T_data yij = ypos[i] - ypos[j];

    CaYY[tid] +=
        ((DPHI_highpass(xij, yij - d, fc, tab_int_x, tab_int_y, Ntab) +
          DPHI_highpass(xij, yij + d, fc, tab_int_x, tab_int_y, Ntab) -
          2 * DPHI_highpass(xij, yij, fc, tab_int_x, tab_int_y, Ntab)) *
         scale * weight);

    tid += blockDim.x * gridDim.x;
  }
}

template __global__ void compute_Ca_element_YY(float *CaYY, int nssp,
                                               float *tab_int_x,
                                               float *tab_int_y, float *xpos,
                                               float *ypos, float d, float fc,
                                               float scale, float weight,
                                               float offset, int Ntab);
template __global__ void compute_Ca_element_YY(
    double *CaYY, int nssp, double *tab_int_x, double *tab_int_y, double *xpos,
    double *ypos, double d, double fc, double scale, double weight,
    double offset, int Ntab);

template <class T_data, typename Fnc, typename Fns>
__device__ void compute_Cerr_element_gen(Fnc const &ptr_cos, Fns const &ptr_sin,
                                         T_data *Cerr, int N, T_data *tab_int_x,
                                         T_data *tab_int_y, T_data *xpos,
                                         T_data *ypos, T_data vdt,
                                         T_data Htheta, T_data L0, T_data fc,
                                         T_data winddir, T_data gsangle,
                                         T_data scale, int Ntab) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < N * N) {
    int i = tid / N;
    int j = tid - i * N;

    T_data xij = xpos[j] - xpos[i];
    T_data yij = ypos[j] - ypos[i];
    T_data xij_vdt = xij - vdt * ptr_cos(winddir);
    T_data yij_vdt = yij - vdt * ptr_sin(winddir);
    T_data xij_ht = xij - Htheta * ptr_cos(gsangle);
    T_data yij_ht = yij - Htheta * ptr_sin(gsangle);
    T_data xij_hvdt = xij_vdt - Htheta * ptr_cos(gsangle);
    T_data yij_hvdt = yij_vdt - Htheta * ptr_sin(gsangle);

    Cerr[tid] +=
        0.5 *
        (DPHI_lowpass(xij_hvdt, yij_hvdt, L0, fc, tab_int_x, tab_int_y, Ntab) -
         DPHI_lowpass(xij_ht, yij_ht, L0, fc, tab_int_x, tab_int_y, Ntab) -
         DPHI_lowpass(xij_vdt, yij_vdt, L0, fc, tab_int_x, tab_int_y, Ntab) +
         DPHI_lowpass(xij, yij, L0, fc, tab_int_x, tab_int_y, Ntab)) *
        scale;
    Cerr[tid] +=
        0.5 *
        (DPHI_lowpass(xij_ht, yij_ht, L0, fc, tab_int_x, tab_int_y, Ntab) -
         DPHI_lowpass(xij, yij, L0, fc, tab_int_x, tab_int_y, Ntab)) *
        scale;
    Cerr[tid] +=
        0.5 *
        (DPHI_lowpass(xij_vdt, yij_vdt, L0, fc, tab_int_x, tab_int_y, Ntab) -
         DPHI_lowpass(xij, yij, L0, fc, tab_int_x, tab_int_y, Ntab)) *
        scale;

    // Cerr[tid] += DPHI_lowpass(xij,yij, L0, fc, tab_int_x, tab_int_y , Ntab);
    tid += blockDim.x * gridDim.x;
  }
}
template <class T_data>
__global__ void compute_Cerr_element(T_data *Cerr, int N, T_data *tab_int_x,
                                     T_data *tab_int_y, T_data *xpos,
                                     T_data *ypos, T_data vdt, T_data Htheta,
                                     T_data L0, T_data fc, T_data winddir,
                                     T_data gsangle, T_data scale, int Ntab,
                                     CarmaDevice *device);

template <>
__global__ void compute_Cerr_element(float *Cerr, int N, float *tab_int_x,
                                     float *tab_int_y, float *xpos, float *ypos,
                                     float vdt, float Htheta, float L0,
                                     float fc, float winddir, float gsangle,
                                     float scale, int Ntab,
                                     CarmaDevice *device) {
  return compute_Cerr_element_gen<float>(cosf, sinf, Cerr, N, tab_int_x,
                                         tab_int_y, xpos, ypos, vdt, Htheta, L0,
                                         fc, winddir, gsangle, scale, Ntab);
}
template <>
__global__ void compute_Cerr_element(double *Cerr, int N, double *tab_int_x,
                                     double *tab_int_y, double *xpos,
                                     double *ypos, double vdt, double Htheta,
                                     double L0, double fc, double winddir,
                                     double gsangle, double scale, int Ntab,
                                     CarmaDevice *device) {
  return compute_Cerr_element_gen<double, double (*)(double),
                                  double (*)(double)>(
      cos, sin, Cerr, N, tab_int_x, tab_int_y, xpos, ypos, vdt, Htheta, L0, fc,
      winddir, gsangle, scale, Ntab);
}

template <class T_data>
__global__ void add_transpose_krnl(T_data *Cerr, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < N * (N + 1) / 2) {
    int ii = N * (N + 1) / 2 - 1 - tid;
    int K = (int)(sqrt(2 * ii + 0.25) - 0.5);
    int i = N - 1 - K;
    int j = tid + i * (i + 1) / 2 - N * i;
    int otid = j + i * N;
    int otidT = i + j * N;
    T_data Cij = Cerr[otid];
    T_data Cji = Cerr[otidT];

    Cerr[otid] = Cij + Cji;
    Cerr[otidT] = Cij + Cji;

    tid += blockDim.x * gridDim.x;
  }
}

template __global__ void add_transpose_krnl<float>(float *Cerr, int N);
template __global__ void add_transpose_krnl<double>(double *Cerr, int N);

template <class T_data>
int compute_Cerr_layer(T_data *Cerr, int N, T_data *tab_int_x,
                       T_data *tab_int_y, T_data *xpos, T_data *ypos,
                       T_data vdt, T_data Htheta, T_data L0, T_data fc,
                       T_data winddir, T_data gsangle, T_data scale, int Ntab,
                       CarmaDevice *device) {
  int nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N * N, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  compute_Cerr_element<<<grid, threads>>>(Cerr, N, tab_int_x, tab_int_y, xpos,
                                          ypos, vdt, Htheta, L0, fc, winddir,
                                          gsangle, scale, Ntab, device);
  carma_check_msg("compute_Cerr_element<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

template int compute_Cerr_layer<float>(float *Cerr, int N, float *tab_int_x,
                                       float *tab_int_y, float *xpos,
                                       float *ypos, float vdt, float Htheta,
                                       float L0, float fc, float winddir,
                                       float gsangle, float scale, int Ntab,
                                       CarmaDevice *device);
template int compute_Cerr_layer<double>(double *Cerr, int N, double *tab_int_x,
                                        double *tab_int_y, double *xpos,
                                        double *ypos, double vdt, double Htheta,
                                        double L0, double fc, double winddir,
                                        double gsangle, double scale, int Ntab,
                                        CarmaDevice *device);

template <class T_data>
int add_transpose(T_data *Cerr, int N, CarmaDevice *device) {
  int nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, N * (N + 1) / 2, nb_blocks, nb_threads);
  dim3 grid2(nb_blocks), threads2(nb_threads);
  add_transpose_krnl<<<grid2, threads2>>>(Cerr, N);
  carma_check_msg("add_transpose<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}
template int add_transpose<float>(float *Cerr, int N, CarmaDevice *device);
template int add_transpose<double>(double *Cerr, int N, CarmaDevice *device);

template <class T_data>
int compute_Ca(T_data *CaXX, T_data *CaYY, int nssp, T_data *tab_int_x,
               T_data *tab_int_y, T_data *xpos, T_data *ypos, T_data offset,
               T_data d, T_data fc, T_data scale, T_data weight, int Ntab,
               CarmaDevice *device) {
  int nb_threads = 0, nb_blocks = 0;
  get_num_blocks_and_threads(device, nssp * nssp, nb_blocks, nb_threads);
  dim3 grid(nb_blocks), threads(nb_threads);

  compute_Ca_element_XX<<<grid, threads>>>(CaXX, nssp, tab_int_x, tab_int_y,
                                           xpos, ypos, d, fc, scale, weight,
                                           offset, Ntab);
  carma_check_msg("compute_Cerr_element_XX<<<>>> execution failed\n");

  compute_Ca_element_YY<<<grid, threads>>>(CaYY, nssp, tab_int_x, tab_int_y,
                                           xpos, ypos, d, fc, scale, weight,
                                           offset, Ntab);
  carma_check_msg("compute_Cerr_element_XX<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

template int compute_Ca(float *CaXX, float *CaYY, int nssp, float *tab_int_x,
                        float *tab_int_y, float *xpos, float *ypos,
                        float offset, float d, float fc, float scale,
                        float weight, int Ntab, CarmaDevice *device);
template int compute_Ca(double *CaXX, double *CaYY, int nssp, double *tab_int_x,
                        double *tab_int_y, double *xpos, double *ypos,
                        double offset, double d, double fc, double scale,
                        double weight, int Ntab, CarmaDevice *device);
