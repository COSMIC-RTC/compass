#include <sutra_atmos.h>
#include <sutra_controller_utils.h>

/*  Tuning parameters of tbulateDPHI kernel*/
#define tabDPHI_thread_x (256)

/*	Tuning parameters of matcov GPU Kernel */
// Thread block size (x, y),
// max #threads per block is 512 for fermi and 1024 for kepler
#define matcov_thread_x (8)
#define matcov_thread_y (8)

//#define CUDA_ERROR_CHECK

//#define CudaSafeCall( err ) __cudaSafeCall( err, __FILE__, __LINE__ )
#define CudaCheckError() __cudaCheckError(__FILE__, __LINE__)
/*
inline void __cudaSafeCall( cudaError err, const char *file, const int line )
{
#ifdef CUDA_ERROR_CHECK
    if ( cudaSuccess != err )
    {
        fprintf( stderr, "cudaSafeCall() failed at %s:%i : %s\n",
                 file, line, cudaGetErrorString( err ) );
        exit( -1 );
    }
#endif

    return;
}
*/
inline void __cudaCheckError(const char *file, const int line) {
#ifdef CUDA_ERROR_CHECK
  cudaError err = cudaGetLastError();
  if (cudaSuccess != err) {
    fprintf(stderr, "cudaCheckError() failed at %s:%i : %s\n", file, line,
            cudaGetErrorString(err));
    exit(-1);
  }

  // More careful checking. However, this will affect performance.
  // Comment away if needed.
  err = cudaDeviceSynchronize();
  if (cudaSuccess != err) {
    fprintf(stderr, "cudaCheckError() with sync failed at %s:%i : %s\n", file,
            line, cudaGetErrorString(err));
    exit(-1);
  }
#endif

  return;
}

//============================================================================================
//================================= AUX FUNCTIONS
//============================================
//============================================================================================
#define VERBOSE 0
void process_err(cudaError_t e, const char *str) {
  if (VERBOSE) printf("%s\n", str);
  if (e != cudaSuccess) {
    printf("*** Error %s: %s \n", str, cudaGetErrorString(e));
    exit(1);
  }
}
//-----------------------------------------------------------------------
double *arr2dAlloc_gpu_gb(long nbLin, long nbCol)
/* DOCUMENT  array = arr2dAlloc(nblin,nbcol)

 Allocates a 2d array (double).
 */
{
  cudaError_t e;
  double *tableau;
  e = cudaMalloc((void **)&tableau, sizeof(double) * nbCol * nbLin);
  process_err(e, "gpu alloc tableau2");
  return tableau;
}

void arr2dFree_gpu_gb(double *tableau)
/* DOCUMENT  arr2dFree(array)

 Free a 2d array (double).
 */
{
  if (tableau) cudaFree(tableau);
}

//============================================================================================
//============================= tabDPHI FUNCTIONS/KERNEL(s)
//==================================
//============================================================================================
__device__ double macdo_x56_gpu_gb(double x, int k)
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
  const double a = 5. / 6.;
  const double x2a = pow(x, (double)2. * a), x22 = x * x / 4.;
  double x2n;  // x^2.a, etc
  double s = 0.0;
  int n;

  const double Ga[11] = {0,
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

  const double Gma[11] = {-3.74878707653729304,     -2.04479295083852408,
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
//------------------------------------------------------------------------------------
__device__ double asymp_macdo_gpu_gb(double x)
/* DOCUMENT asymp_macdo_gpu_gb(x)

 Computes a term involved in the computation of the phase struct
 function with a finite outer scale according to the Von-Karman
 model. The term involves the MacDonald function (modified bessel
 function of second kind) K_{5/6}(x), and the algorithm uses the
 asymptotic form for x ~ infinity.
 Warnings :
 - This function makes a doubleing point interrupt for x=0
 and should not be used in this case.
 - Works only for x>0.

 SEE ALSO:
 */
{
  // k2 is the value for
  // gamma_R(5./6)*2^(-1./6)
  const double k2 = 1.00563491799858928388289314170833;
  const double k3 = 1.25331413731550012081;   //  sqrt(pi/2)
  const double a1 = 0.22222222222222222222;   //  2/9
  const double a2 = -0.08641975308641974829;  //  -7/89
  const double a3 = 0.08001828989483310284;   // 175/2187
  double res;
  double x_1;

  x_1 = 1. / x;
  res = k2 - k3 * exp(-x) * pow(x, (double)(1 / 3.)) *
                 (1.0 + x_1 * (a1 + x_1 * (a2 + x_1 * a3)));
  return res;
}
//------------------------------------------------------------------------------------
__device__ double rodconan_gpu_gb(double r, double L0, int k)
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
  const double pi = 3.1415926535897932384626433;
  double res = 0;

  // k1 is the value of :
  // 2*gamma_R(11./6)*2^(-5./6)*pi^(-8./3)*(24*gamma_R(6./5)/5.)^(5./6);
  const double k1 = 0.1716613621245709486;
  const double dprf0 = (2 * pi / L0) * r;
  // k2 is the value for gamma_R(5./6)*2^(-1./6),
  // but is now unused
  // k2 = 1.0056349179985892838;

  // Xlim = 0.75*2*pi;   // = 4.71239
  if (dprf0 > 4.71239)
    res = asymp_macdo_gpu_gb(dprf0);
  else
    res = -macdo_x56_gpu_gb(dprf0, k);

  res *= k1 * pow(L0, (double)5. / 3);

  return res;
}

__global__ void tabulateDPHI_gpu_gb_kernel(double *tabDPHI_d, double *L0diff_d,
                                           long Nl0, long Ndphi,
                                           double convert) {
  const int tx = threadIdx.x;
  const int ty = blockIdx.x;

  const int tid = ty * blockDim.x + tx;
  int l = tid / Ndphi;
  int j = tid % Ndphi;

  if (tid >= (Nl0 * Ndphi)) return;

  tabDPHI_d[tid] = rodconan_gpu_gb((double)j / convert, L0diff_d[l], 10);

  // double* mytabDPHI = tabDPHI_d + (l * Ndphi);
  //
  // int j, k;
  //#pragma unroll
  // for(k = 0; k < (Ndphi/tabDPHI_thread_x); k++)
  //{
  //	j = k * tabDPHI_thread_x + tx;
  //	mytabDPHI[j] = rodconan_gpu_gb(rr_d[j], L0diff_d[l], 10);
  //}
  //
  // k = (Ndphi/tabDPHI_thread_x);
  // if(tx < (Ndphi%tabDPHI_thread_x) )
  //{
  //	j = k * tabDPHI_thread_x + tx;
  //	mytabDPHI[j] = rodconan_gpu_gb(rr_d[j], L0diff_d[l], 10);
  //}
}

__device__ double Ij0t83_gb(double x, double *tab_x, double *tab_y, long npts) {
  if (x <= exp(-3.0))
    return (double)(0.75 * pow((double)x, 1 / 3.) * (1 - x * x / 112.));
  else {
    double dt = 14. / (npts - 1);  // 14 = tmax - tmin
    double convert = (log(x) + 4.) / dt;
    long i0 = (long)convert;
    long i1 = i0 + 1;
    /*
    long i1 = 0;
    while(x > tab_x[i1])
        i1++;
    long i0 = i1 - 1;
    */
    // long i0 = (long) ((logf(x) + 4.0f) * convert);
    // long i1 = i0 + 1;
    // return  ((x - (double)i0 / convert) * tab_y[indexL0 * npts + i1]
    //    + ((double)i1 / convert - x) * tab_y[indexL0 * npts + i0]);
    return (((x - tab_x[i0]) * tab_y[/*indexL0 * npts +*/ i1] +
             (tab_x[i1] - x) * tab_y[/*indexL0 * npts +*/ i0])) /
           (tab_x[i1] - tab_x[i0]);
  }
}

__global__ void tabulateDPHI_lowpass_kernel(double *tabDPHI_d,
                                            double *tab_int_x,
                                            double *tab_int_y, double *dx,
                                            double *L0diff_d, long Nl0,
                                            long Ndphi, double convert,
                                            double convert_int, long npts) {
  const double pi = 3.1415926535897932384626433;
  const int tx = threadIdx.x;
  const int ty = blockIdx.x;

  const int tid = ty * blockDim.x + tx;
  int dm = tid / (Ndphi * Nl0);
  int pos = tid - dm * (Ndphi * Nl0);
  int l = pos / Ndphi;
  // int j = tid % Ndphi;
  int j = pos - l * Ndphi;
  double r = (double)j / convert;
  // long indexL0 = 0;

  if (tid >= (Nl0 * Ndphi)) return;

  // tabDPHI_d[tid] = Ij0t83_gb((double)(r*(pi/dx)), indexL0, tab_int_x,
  // tab_int_y, convert_int,npts);

  tabDPHI_d[tid] =
      pow(r, (double)(5. / 3.)) *
      Ij0t83_gb((double)(r * (pi / dx[dm])), tab_int_x, tab_int_y, npts) *
      (double)((2 * pow((2 * pi), (double)(8 / 3.)) * 0.0228956));
}

//------------------------------------------------------------------------------------
__device__ double DPHI_gpu_gb(double x, double y, double L0)
/* DOCUMENT dphi = DPHI(x,y,indexL0,rr,tabDPHI,convert) * r0^(-5./3)
 <x> & <y>         :  separation between apertures
 <indexL0>         :  index for the L0 taken into account
 <rr>              :  array of distance between apertures
 <tabDPHI>         :  array of precomputed DPHI
 <convert>         :  relation between the index on tabDPHI and (x,y)

 Computes the phase structure function for a separation (x,y).
 The r0 is not taken into account : the final result of DPHI(x,y,L0)
 has to be scaled with r0^-5/3, with r0 expressed in meters, to get
 the right value.

 SEE ALSO:
 */
{
  double r = sqrt(x * x + y * y);

  return rodconan_gpu_gb(r, L0, 10);
  /*
  long i0 = (long) (r * convert);
  long i1 = i0 + 1;

  return ((r - (double)i0 / convert) * tabDPHI[indexL0 * Ndphi + i1]
    + ((double)i1 / convert - r) * tabDPHI[indexL0 * Ndphi + i0]);
    */
}

//============================================================================================
//============================= SUBAP POSITION KERNELS/FUNCTIONS
//=============================
//============================================================================================
__global__ void subposition_gpu_gb_kernel(
    long Nw, long Nx, long *Nsubap, long Nlayer, double *alphaX, double *alphaY,
    double *h, double *GsAlt, long *Nssp, double *diamPup, double *thetaML,
    long *ioff, double *X, double *Y, double *XPup, double *YPup, double *u,
    double *v) {
  const int tx = threadIdx.x;
  const int ty = blockIdx.x;

  const int tid = ty * blockDim.x + tx;
  long i;      // subaperture i
  long n = 0;  // WFS n
  long l;
  const double rad = 3.14159265358979323846 / 180.;

  if (tid >= (Nx * Nlayer)) return;

  l = tid / Nx;

  const int pos = tid - l * Nx;
  long Nsubapx = Nsubap[0];

  while (pos >= Nsubapx) {
    n++;
    Nsubapx += Nsubap[n];
  }
  Nsubapx -= Nsubap[n];

  i = pos - Nsubapx;

  // tid = n + i * Nw + l * Nw * Nsubap

  const double dX = alphaX[n] * h[l];
  const double dY = alphaY[n] * h[l];

  const double rr = 1. - h[l] * GsAlt[n];

  const long nssp = Nssp[n];

  // magnification factor
  const double G = diamPup[n] / (double)(nssp);

  // rotation angle
  const double th = thetaML[n] * rad;

  // taking magnification factor into account
  const double xtp = X[ioff[n] + i] * G;
  const double ytp = Y[ioff[n] + i] * G;

  // taking rotation into account
  double uu = xtp * cos(th) - ytp * sin(th);
  double vv = xtp * sin(th) + ytp * cos(th);

  // taking pupil offset into account
  uu += XPup[n];
  vv += YPup[n];

  // Projection onto  the layer

  u[tid] = uu * rr + dX;
  v[tid] = vv * rr + dY;
}

//============================================================================================
//============================= MATCOV ELEMENTARY FUNCTIONS
//==================================
//============================================================================================
__device__ double cov_XX_gpu_gb(double du, double dv, double ac, double ad,
                                double bc, double bd, double *tabDPHI,
                                double L0, double convert, int Ndphi)
/* DOCUMENT
  Compute the XX-covariance with the distance sqrt(du2+dv2). DPHI is precomputed
  on tabDPHI.
*/
{
  return -DPHI_gpu_gb(du + ac, dv, L0) + DPHI_gpu_gb(du + ad, dv, L0) +
         DPHI_gpu_gb(du + bc, dv, L0) - DPHI_gpu_gb(du + bd, dv, L0);
}

//------------------------------------------------------------------------------------
__device__ double cov_YY_gpu_gb(double du, double dv, double ac, double ad,
                                double bc, double bd, double *tabDPHI,
                                double L0, double convert, int Ndphi)
/* DOCUMENT
   Compute the YY-covariance with the distance sqrt(du2+dv2). DPHI is
   precomputed on tabDPHI.
 */
{
  return -DPHI_gpu_gb(du, dv + ac, L0) + DPHI_gpu_gb(du, dv + ad, L0) +
         DPHI_gpu_gb(du, dv + bc, L0) - DPHI_gpu_gb(du, dv + bd, L0);
}

//------------------------------------------------------------------------------------
__device__ double cov_XY_gpu_gb(double du, double dv, double s0,
                                double *tabDPHI, double L0, double convert,
                                int Ndphi)
/* DOCUMENT
   Compute the XY-covariance with the distance sqrt(du2+dv2). DPHI is
   precomputed on tabDPHI.
 */
{
  return -DPHI_gpu_gb(du + s0, dv - s0, L0) +
         DPHI_gpu_gb(du + s0, dv + s0, L0) + DPHI_gpu_gb(du - s0, dv - s0, L0) -
         DPHI_gpu_gb(du - s0, dv + s0, L0);
}

//============================================================================================
//============================= CPHIM ELEMENTARY FUNCTIONS
//==================================
//============================================================================================
__device__ double DPHI_highpass_gb(double r, double fc, double *tab_x,
                                   double *tab_y, long npts) {
  const double pi = 3.1415926535897932384626433;
  return pow(r, 5 / 3.) *
         (1.1183343328701949 - Ij0t83_gb(2 * pi * fc * r, tab_x, tab_y, npts)) *
         pow(2 * pi, 8 / 3.) * 2 * 0.0228956;
}
__device__ double DPHI_lowpass_gb(double x, double y, double L0, double fc,
                                  double *tab_int_x, double *tab_int_y,
                                  long npts) {
  /*
  double r = sqrt(x * x + y * y);
  const double pi = 3.1415926535897932384626433;
  int npts = (int)1/pas;
  double du = 2*pi*fc*r/npts;
  double dphi = 0;
  double u = 0;
  //for(double u=pas ; u <= 2*pi*fc*r ; u+=pas ){
  for (int i=1 ; i < npts ; i++){
        u += du;
        dphi += pow(u*u + pow(2*pi*r/L0,2),-11/6.) * u * (1- j0(u)) * du;
  }

  return 2*pow((2*pi),8/3.)*0.0228956*pow(r,5/3.) * dphi;
  */
  double r = sqrt(x * x + y * y);

  return rodconan_gpu_gb(r, L0, 10) -
         DPHI_highpass_gb(r, fc, tab_int_x, tab_int_y, npts);
}
__device__ double cphim_XX(double du, double dv, double posx, double posy,
                           double xref, double yref, double s2, double L0,
                           double fc, long npts, double *tab_int_x,
                           double *tab_int_y)
/* DOCUMENT
  Compute the XX-covariance with the distance sqrt(du2+dv2). DPHI is precomputed
  on tabDPHI.
*/
{
  return -DPHI_lowpass_gb(du - 2 * s2, dv - s2, L0, fc, tab_int_x, tab_int_y,
                          npts) +
         DPHI_lowpass_gb(du, dv - s2, L0, fc, tab_int_x, tab_int_y, npts) +
         DPHI_lowpass_gb(posx + 2 * s2 - xref, posy + s2 - yref, L0, fc,
                         tab_int_x, tab_int_y, npts) -
         DPHI_lowpass_gb(posx - xref, posy + s2 - yref, L0, fc, tab_int_x,
                         tab_int_y, npts);

  /*
                return -DPHI_gpu_gb(du - 2*s2, dv - s2, L0)
                                + DPHI_gpu_gb(du, dv - s2, L0)
                                + DPHI_gpu_gb(posx+ 2*s2, posy + s2, L0)
                                - DPHI_gpu_gb(posx, posy+s2, L0);*/
}

//------------------------------------------------------------------------------------
__device__ double cphim_YY(double du, double dv, double posx, double posy,
                           double xref, double yref, double s2, double L0,
                           double fc, long npts, double *tab_int_x,
                           double *tab_int_y)
/* DOCUMENT
   Compute the YY-covariance with the distance sqrt(du2+dv2). DPHI is
   precomputed on tabDPHI.
 */
{
  return -DPHI_lowpass_gb(du - s2, dv - 2 * s2, L0, fc, tab_int_x, tab_int_y,
                          npts) +
         DPHI_lowpass_gb(du - s2, dv, L0, fc, tab_int_x, tab_int_y, npts) +
         DPHI_lowpass_gb(posx + s2 - xref, posy + 2 * s2 - yref, L0, fc,
                         tab_int_x, tab_int_y, npts) -
         DPHI_lowpass_gb(posx + s2 - xref, posy - yref, L0, fc, tab_int_x,
                         tab_int_y, npts);

  /*
    return  -DPHI_gpu_gb(du-s2, dv - 2*s2, L0)
      + DPHI_gpu_gb(du-s2, dv, L0)
      + DPHI_gpu_gb(posx+s2, posy + 2*s2, L0)
      - DPHI_gpu_gb(posx+s2, posy, L0);*/
}
//============================================================================================
//============================= MATCOV 3 FUNCTIONS/KERNEL
//====================================
//============================================================================================
__device__ double compute_element_3(int ipos, int jpos, double convert,
                                    double *sspSizeL, long *Nssp, double *u,
                                    double *v, double pasDPHI, double *tabDPHI,
                                    long *indexL0, double *cn2, int Ndphi,
                                    int Nw, int Nlayer, int Nsubap,
                                    int type_mat, double teldiam) {
  /* *** Covariance matrix per-element generation ***
   *   Arguments
   *   =========
   *	ipos:		Integer: global x-coordinate of the element w.r.t. the
   *entire matrix jpos:		Integer: global y-coordinate of the element
   *w.r.t. the entire matrix
   */

  const double lambda2 = 0.00026942094446267851;
  const int nslps = Nsubap * 2;

  // WFS m
  int m = ipos / nslps;  // tab_wfs[ipos];
  if (type_mat == 3) m = Nw - 1;
  // WFS n
  int n = jpos / nslps;  // tab_wfs[jpos];
  if (type_mat == 2) n = Nw - 1;

  // subap i
  int i = ipos % (nslps / 2);  // tab_subap[ipos];
  // subap j
  int j = jpos % (nslps / 2);  // tab_subap[jpos];;

  // xy i
  int xy_i = (ipos / (nslps / 2)) % 2;  // tab_xy[ipos];
  // xy j
  int xy_j = (jpos / (nslps / 2)) % 2;  // tab_xy[jpos];

  const double sspSizem = teldiam / Nssp[m];
  const double sspSizen = teldiam / Nssp[n];

  const double kk = lambda2 / (sspSizem * sspSizen);

  int type = xy_i * 2 + xy_j;

  // Layer l
  double covar = 0.0;
#pragma unroll
  for (int l = 0; l < Nlayer; l++) {
    const double sspSizeml = sspSizeL[m * Nlayer + l];
    const double sspSizenl = sspSizeL[n * Nlayer + l];
    // test if the altitude layers is not higher than the LGS altitude
    if ((sspSizeml > 0) && (sspSizenl > 0)) {
      const int pos1 = m + i * Nw + l * Nw * Nsubap;
      const int pos2 = n + j * Nw + l * Nw * Nsubap;
      const double du = u[pos1] - u[pos2];
      const double dv = v[pos1] - v[pos2];

      const double s1 = sspSizeml * 0.5;
      const double s2 = sspSizenl * 0.5;

      const double ac = s1 - s2;
      const double ad = s1 + s2;
      const double bc = -ad;  // initially -s1-s2;
      const double bd = -ac;  // initially -s1+s2;

      if (type == 0)
        covar += 0.5 * pasDPHI *
                 cov_XX_gpu_gb(du, dv, ac, ad, bc, bd, tabDPHI, indexL0[l],
                               convert, Ndphi) *
                 kk * cn2[l];
      else if (type == 3)
        covar += 0.5 * pasDPHI *
                 cov_YY_gpu_gb(du, dv, ac, ad, bc, bd, tabDPHI, indexL0[l],
                               convert, Ndphi) *
                 kk * cn2[l];
      else {  // if ((type == 1) || (type == 2))
        const double s0 =
            sqrt(s1 * s1 + s2 * s2);  // half size of the subaperture equivalent
                                      // to a convolution by s1 and s2
        const double dd =
            (s1 > s2) ? 1. - s2 / s1 : 1. - s1 / s2;  // Nono's style ....
        covar +=
            0.25 * pasDPHI *
            cov_XY_gpu_gb(du, dv, s0, tabDPHI, indexL0[l], convert, Ndphi) *
            kk * cn2[l] * (1. - dd * dd);
      }
    }
  }
  return (double)covar;
}

__global__ void matcov_kernel_3(char uplo, char copy, double *data, int nrows,
                                int ncols, int xoffset, int yoffset, int lda,
                                double convert, double *sspSizeL, long *Nssp,
                                double *u, double *v, double pasDPHI,
                                double *tabDPHI, long *indexL0, double *cn2,
                                int Ndphi, int Nw, int Nlayer, int Nsubap,
                                int type_mat, double teldiam) {
  /* *** covariance matrix generation kernel ***
   *	The kernel generates the element values in a given matrix/submatrix
   *   The generation function can be any function, as long as each element
   *   can be computed both individually and independently
   *
   *	see argument description in the kernel driver
   */

  // local thread coordinates w.r.t. thread block
  const int tx_ = threadIdx.x;
  const int ty_ = threadIdx.y;

  // local thread block coordinates w.r.t. kernel grid
  const int bx_ = blockIdx.x;
  const int by_ = blockIdx.y;

  // local coordinates of the element w.r.t. submatrix
  int lx = bx_ * blockDim.x + tx_;
  int ly = by_ * blockDim.y + ty_;

  // global coordinates of the elemnt w.r.t. the entire matrix
  int gx = lx + xoffset;
  int gy = ly + yoffset;

  // out-of-bound threads should terminate
  if ((lx >= nrows) || (ly >= ncols)) return;

  double value;
  if (uplo == 'l') {
    if (gy <= gx) {
      value = compute_element_3(gx, gy, convert, sspSizeL, Nssp, u, v, pasDPHI,
                                tabDPHI, indexL0, cn2, Ndphi, Nw, Nlayer,
                                Nsubap, type_mat, teldiam);
      data[ly * lda + lx] = value;
      if (copy == 'c') data[lx * lda + ly] = value;
    }
  } else if (uplo == 'u') {  // upper
    if (gx <= gy) {
      value = compute_element_3(gx, gy, convert, sspSizeL, Nssp, u, v, pasDPHI,
                                tabDPHI, indexL0, cn2, Ndphi, Nw, Nlayer,
                                Nsubap, type_mat, teldiam);
      data[ly * lda + lx] = value;
      if (copy == 'c') data[lx * lda + ly] = value;
    }
  } else {  // uplo = 'f' full generation
    value = compute_element_3(gx, gy, convert, sspSizeL, Nssp, u, v, pasDPHI,
                              tabDPHI, indexL0, cn2, Ndphi, Nw, Nlayer, Nsubap,
                              type_mat, teldiam);
    data[ly * lda + lx] = value;
  }

  // if ((type_mat == 3) || (gx <= gy))
  //{
  // call the generation function
  // data[0] = compute_element_3(gx, gy, tab_wfs, tab_subap,
  // tab_xy,convert,sspSizeL,Nssp,u,v,pasDPHI,tabDPHI,
  //		      indexL0,cn2,Ndphi,Nw,Nlayer,Nsubap,type_mat,teldiam);
  // printf("gx = %d, gy = %d ----- %.2f \n", gx, gy, data[0]);
  //}
}

//============================================================================================
//============================= MATCOV TS FUNCTIONS/KERNEL
//===================================
//============================================================================================
__device__ double compute_element_ts_(int ipos, int jpos, double convert,
                                      double *X, double *Y, long *Nssp,
                                      double pasDPHI, double *tabDPHI,
                                      long *indexL0, double *cn2, int Ndphi,
                                      int Nw, int Nlayer, int Nsubap,
                                      double teldiam) {
  /* *** Covariance matrix per-element generation ***
   *   Arguments
   *   =========
   *	ipos:		Integer: global x-coordinate of the element w.r.t. the
   *entire matrix jpos:		Integer: global y-coordinate of the element
   *w.r.t. the entire matrix
   */

  // for now return a dummy value

  double lambda2 = 0.00026942094446267851;
  // WFS Nw-1
  // subap i
  int i = ipos < Nsubap ? ipos : ipos - Nsubap;
  // subap j
  int j = jpos < Nsubap ? jpos : jpos - Nsubap;
  // xy i
  int xy_i = ipos < Nsubap ? 0 : 1;
  // xy j
  int xy_j = jpos < Nsubap ? 0 : 1;

  double sspSize = teldiam / Nssp[Nw - 1];

  double kk = lambda2 / (sspSize * sspSize);

  int type = xy_i * 2 + xy_j;

  double s = sspSize * 0.5;

  double ac = 0.0;
  double ad = 2.0 * s;
  double bc = -ad;
  double bd = 0.0;

  double du = X[(Nsubap * (Nw - 1) + i)] - X[(Nsubap * (Nw - 1) + j)];
  double dv = Y[(Nsubap * (Nw - 1) + i)] - Y[(Nsubap * (Nw - 1) + j)];

  // if(ipos < 10)printf("ipos = %d - %d\n", ipos, (Nsubap*(Nw-1)+i));
  // if(jpos < 10)printf("jpos = %d - %d\n", jpos, (Nsubap*(Nw-1)+j));

  // const double du = X[0] - X[1];
  // const double dv = Y[0] - Y[1];

  // Layer l
  double covar = 0.0;
#pragma unroll
  for (int l = 0; l < Nlayer; l++) {
    // test if the altitude layers is not higher than the LGS altitude
    if (sspSize > 0) {
      if (type == 0)
        covar += 0.5 * pasDPHI *
                 cov_XX_gpu_gb(du, dv, ac, ad, bc, bd, tabDPHI, indexL0[l],
                               convert, Ndphi) *
                 kk * cn2[l];
      else if (type == 3)
        covar += 0.5 * pasDPHI *
                 cov_YY_gpu_gb(du, dv, ac, ad, bc, bd, tabDPHI, indexL0[l],
                               convert, Ndphi) *
                 kk * cn2[l];
      else {
        double s0 = 1.41421 * s;  // half size of the subaperture equivalent to
                                  // a convolution by s1 and s2
        double dd = 0;
        covar +=
            0.25 * pasDPHI *
            cov_XY_gpu_gb(du, dv, s0, tabDPHI, indexL0[l], convert, Ndphi) *
            kk * cn2[l] * (1. - dd * dd);
      }
    }
  }
  return (double)covar;
}
//--------------------------------------------------------------------------------------------
__global__ void matcov_ts_kernel(double *data, int nrows, int ncols,
                                 int xoffset, int yoffset, int lda,
                                 double convert, double *X, double *Y,
                                 long *Nssp, double pasDPHI, double *tabDPHI,
                                 long *indexL0, double *cn2, int Ndphi, int Nw,
                                 int Nlayer, int Nsubap, double teldiam) {
  /* *** covariance matrix generation kernel ***
   *	The kernel generates the element values in a given matrix/submatrix
   *   The generation function can be any function, as long as each element
   *   can be computed both individually and independently
   *
   *	see argument description in the kernel driver
   */

  // local thread coordinates w.r.t. thread block
  const int tx_ = threadIdx.x;
  const int ty_ = threadIdx.y;

  // local thread block coordinates w.r.t. kernel grid
  const int bx_ = blockIdx.x;
  const int by_ = blockIdx.y;

  // local coordinates of the element w.r.t. submatrix
  int lx = bx_ * blockDim.x + tx_;
  int ly = by_ * blockDim.y + ty_;

  // global coordinates of the elemnt w.r.t. the entire matrix
  int gx = lx + xoffset;
  int gy = ly + yoffset;

  // out-of-bound threads should terminate
  if ((lx >= nrows) || (ly >= ncols)) return;

  // Advance the data pointer accordingly
  data += ly * lda + lx;

  // call the generation function
  data[0] =
      compute_element_ts_(gx, gy, convert, X, Y, Nssp, pasDPHI, tabDPHI,
                          indexL0, cn2, Ndphi, Nw, Nlayer, Nsubap, teldiam);
  // printf("gx = %d, gy = %d ----- %.2f \n", gx, gy, data[0]);
}

//============================================================================================
//============================= MATCOV TS
//===================================================
//============================================================================================
/*
void matts_gpu_gb(double* data, int nrows, int ncols, int xoffset, int yoffset,
int lda, struct tomo_struct tomo, struct gtomo_struct *tomo_gpu)
{
        /* *** matcov gpu kernel driver ***
        *  Arguments
        *  ==========
        *  data		double pointer: A pointer to the matrix/submatrix to be
generated. It
        *  			should always point to the first element in a
matrix/submatrix
        *
        *  nrows	integer: The number of rows of the matrix/submatrix to
be generated
        *
        *  ncols	integer: The number of columns of the matrix/submatrix
to be generated
        *
        *  xoffset	integer: The x-offset of the submatrix, must be zero if
the entire matrix
        *			is generated. Its the x-coordinate of the first
element in the matrix/submatrix
        *
        *  yoffset  integer: The y-offset of the submatrix, must be zero if the
entire matrix
        *			is generated. Its the y-coordinate of the first
element in the matrix/submatrix
        *
        *  lda		integer: The leading dimension of the matrix/submatrix
        */
/*
  const long Nw = tomo.Nw;
  const double crmax = tomo.rmax;
  const double pasDPHI = 1./tomo.pasDPHI; //inverse du pas de rr
  const long Ndphi = floor(crmax*pasDPHI)+1;
  const double convert = (double)(Ndphi-1)/(crmax+1./pasDPHI);


  int nbx = nrows / matcov_thread_x + (nrows%matcov_thread_x != 0);
  int nby = ncols / matcov_thread_y + (ncols%matcov_thread_y != 0);

  dim3 dimBlock(matcov_thread_x, matcov_thread_y);
  dim3 dimGrid(nbx, nby);
  const long Nsubap = tomo.Nsubap[Nw-1];

  matcov_ts_kernel<<<dimGrid, dimBlock, 0, tomo_gpu->matcov_stream>>>(data,
nrows, ncols, xoffset, yoffset, lda,
                                           convert,tomo_gpu->X_d,tomo_gpu->Y_d,tomo_gpu->Nssp_d,
                                           pasDPHI,tomo_gpu->tabDPHI_d,tomo_gpu->indexL0_d,tomo_gpu->cn2_d,
                                           Ndphi,tomo.Nw,atmos->nscreens,Nsubap,tomo.DiamTel);
  //CudaCheckError();
}
*/
//============================================================================================
//============================= MATCOV COPY KERNEL
//===========================================
//============================================================================================
__global__ void matcov_kernel_copy(double *data, int nrows, int ncols,
                                   int xoffset, int yoffset, int lda) {
  /* *** covariance matrix generation kernel ***
   *	The kernel generates the element values in a given matrix/submatrix
   *   The generation function can be any function, as long as each element
   *   can be computed both individually and independently
   *
   *	see argument description in the kernel driver
   */

  // local thread coordinates w.r.t. thread block
  const int tx_ = threadIdx.x;
  const int ty_ = threadIdx.y;

  // local thread block coordinates w.r.t. kernel grid
  const int bx_ = blockIdx.x;
  const int by_ = blockIdx.y;

  // local coordinates of the element w.r.t. submatrix
  int lx = bx_ * blockDim.x + tx_;
  int ly = by_ * blockDim.y + ty_;

  // global coordinates of the elemnt w.r.t. the entire matrix
  int gx = lx + xoffset;
  int gy = ly + yoffset;

  // out-of-bound threads should terminate
  if ((lx >= nrows) || (ly >= ncols)) return;

  // Advance the data pointer accordingly
  // data += ly * lda + lx;

  if (gx > gy) {
    // call the generation function
    data[ly * lda + lx] = data[ly + lx * lda];
    // printf("gx = %d, gy = %d ----- %.2f \n", gx, gy, data[0]);
  }
}

//============================================================================================
//============================= MATCOV 1
//=====================================================
//============================================================================================
//************************** OBSOLETE - REMOVED
//********************************************//

//============================================================================================
//============================= MATCOV 2
//=====================================================
//============================================================================================
//************************** OBSOLETE - REMOVED
//********************************************//

//============================================================================================
//=============================== TOMO INIT/FIN FUNCTIONS
//====================================
//============================================================================================
void init_tomo_gpu_gb(struct gtomo_struct *tomo_gpu, sutra_atmos *atmos,
                      sutra_sensors *sensors, double diamTel, double cobs) {
  cudaError_t e;

  tomo_gpu->DiamTel = diamTel;
  tomo_gpu->obs = cobs;
  tomo_gpu->Nw = sensors->nsensors();  // Adding TS for debug

  tomo_gpu->lgs_cst = 0.;
  tomo_gpu->spot_width = 1.;
  tomo_gpu->lgs_depth = 10000.;
  tomo_gpu->lgs_alt = 90000.;
  tomo_gpu->nlgs = 0;
  tomo_gpu->pasDPHI = 0.0001;

  tomo_gpu->Nx = 0;
  tomo_gpu->Nssp = (long *)malloc(tomo_gpu->Nw * sizeof(long));
  tomo_gpu->Nsubap = (long *)malloc(tomo_gpu->Nw * sizeof(long));
  tomo_gpu->diamPup = (double *)malloc(tomo_gpu->Nw * sizeof(double));
  tomo_gpu->XPup = (double *)malloc(tomo_gpu->Nw * sizeof(double));
  tomo_gpu->YPup = (double *)malloc(tomo_gpu->Nw * sizeof(double));
  tomo_gpu->thetaML = (double *)malloc(tomo_gpu->Nw * sizeof(double));
  tomo_gpu->GsAlt = (double *)malloc(tomo_gpu->Nw * sizeof(double));

  for (int i = 0; i < tomo_gpu->Nw; i++) {
    tomo_gpu->Nssp[i] = sensors->d_wfs[i]->nxsub;
    tomo_gpu->Nsubap[i] = sensors->d_wfs[i]->nvalid;
    tomo_gpu->diamPup[i] = (double)tomo_gpu->Nssp[i];
    tomo_gpu->XPup[i] = 0.;
    tomo_gpu->YPup[i] = 0.;
    tomo_gpu->thetaML[i] = 0.;
    if (sensors->d_wfs[i]->d_gs->lgs) {
      tomo_gpu->nlgs += 1;
      tomo_gpu->GsAlt[i] = 1.0 / tomo_gpu->lgs_alt;
    } else
      tomo_gpu->GsAlt[i] = 0.0;
    tomo_gpu->Nx += sensors->d_wfs[i]->nvalid;
  }

  e = cudaMalloc((void **)&(tomo_gpu->indexL0_d),
                 atmos->nscreens * sizeof(long));
  process_err(e, "alloc gpu indexL0_d");

  e = cudaMalloc((void **)&(tomo_gpu->u_d),
                 atmos->nscreens * tomo_gpu->Nx * sizeof(double));
  process_err(e, "alloc gpu u_d");

  e = cudaMalloc((void **)&(tomo_gpu->v_d),
                 atmos->nscreens * tomo_gpu->Nx * sizeof(double));
  process_err(e, "alloc gpu v_d");

  e = cudaMalloc((void **)&(tomo_gpu->sspSizeL_d),
                 tomo_gpu->Nw * atmos->nscreens * sizeof(double));
  process_err(e, "alloc gpu sspSizeL_d");

  e = cudaMalloc((void **)&(tomo_gpu->cn2_d), atmos->nscreens * sizeof(double));
  process_err(e, "alloc gpu cn2_d");

  e = cudaMalloc((void **)&(tomo_gpu->h_d), atmos->nscreens * sizeof(double));
  process_err(e, "alloc gpu h_d");

  e = cudaMalloc((void **)&(tomo_gpu->Nssp_d), tomo_gpu->Nw * sizeof(long));
  process_err(e, "alloc gpu Nssp_d");

  e = cudaMalloc((void **)&(tomo_gpu->Nsubap_d), tomo_gpu->Nw * sizeof(long));
  process_err(e, "alloc gpu Nsubap_d");

  e = cudaMalloc((void **)&(tomo_gpu->ioff_d), tomo_gpu->Nw * sizeof(long));
  process_err(e, "alloc gpu ioff_d");

  e = cudaMalloc((void **)&(tomo_gpu->alphaX_d), tomo_gpu->Nw * sizeof(double));
  process_err(e, "alloc gpu alphaX_d");

  e = cudaMalloc((void **)&(tomo_gpu->alphaY_d), tomo_gpu->Nw * sizeof(double));
  process_err(e, "alloc gpu alphaY_d");

  e = cudaMalloc((void **)&(tomo_gpu->GsAlt_d), tomo_gpu->Nw * sizeof(double));
  process_err(e, "alloc gpu GsAlt_d");

  e = cudaMalloc((void **)&(tomo_gpu->diamPup_d),
                 tomo_gpu->Nw * sizeof(double));
  process_err(e, "alloc gpu diamPup_d");

  e = cudaMalloc((void **)&(tomo_gpu->thetaML_d),
                 tomo_gpu->Nw * sizeof(double));
  process_err(e, "alloc gpu thetaML_d");

  e = cudaMalloc((void **)&(tomo_gpu->X_d), tomo_gpu->Nx * sizeof(double));
  process_err(e, "alloc gpu X_d");

  e = cudaMalloc((void **)&(tomo_gpu->Y_d), tomo_gpu->Nx * sizeof(double));
  process_err(e, "alloc gpu Y_d");

  e = cudaMalloc((void **)&(tomo_gpu->XPup_d), tomo_gpu->Nw * sizeof(double));
  process_err(e, "alloc gpu XPup_d");

  e = cudaMalloc((void **)&(tomo_gpu->YPup_d), tomo_gpu->Nw * sizeof(double));
  process_err(e, "alloc gpu YPup_d");

  tomo_gpu->L0diff_d = NULL;
  tomo_gpu->tabDPHI_d = NULL;

  e = cudaStreamCreate(&(tomo_gpu->matcov_stream));
  process_err(e, "create matcov stream");
}

void free_tomo_gpu_gb(struct gtomo_struct *tomo_gpu) {
  cudaError_t e;

  if ((tomo_gpu->u_d)) e = cudaFree(tomo_gpu->u_d);
  process_err(e, "free gpu u_d");

  if (tomo_gpu->v_d) e = cudaFree(tomo_gpu->v_d);
  process_err(e, "free gpu v_d");

  if (tomo_gpu->sspSizeL_d) e = cudaFree(tomo_gpu->sspSizeL_d);
  process_err(e, "free gpu sspSizeL_d");

  if (tomo_gpu->cn2_d) e = cudaFree(tomo_gpu->cn2_d);
  process_err(e, "free gpu cn2_d");

  if (tomo_gpu->h_d) e = cudaFree(tomo_gpu->h_d);
  process_err(e, "free gpu h_d");

  if (tomo_gpu->Nsubap_d) e = cudaFree(tomo_gpu->Nsubap_d);
  process_err(e, "free gpu Nsubap_d");

  if (tomo_gpu->indexL0_d) e = cudaFree(tomo_gpu->indexL0_d);
  process_err(e, "free gpu indexL0_d");

  if (tomo_gpu->Nssp_d) e = cudaFree(tomo_gpu->Nssp_d);
  process_err(e, "free gpu Nssp_d");

  if (tomo_gpu->ioff_d) e = cudaFree(tomo_gpu->ioff_d);
  process_err(e, "free gpu ioff_d");

  if (tomo_gpu->alphaX_d) e = cudaFree(tomo_gpu->alphaX_d);
  process_err(e, "free gpu alphaX_d");

  if (tomo_gpu->alphaY_d) e = cudaFree(tomo_gpu->alphaY_d);
  process_err(e, "free gpu alphaY_d");

  if (tomo_gpu->GsAlt_d) e = cudaFree(tomo_gpu->GsAlt_d);
  process_err(e, "free gpu GsAlt_d");

  if (tomo_gpu->diamPup_d) e = cudaFree(tomo_gpu->diamPup_d);
  process_err(e, "free gpu diamPup_d");

  if (tomo_gpu->thetaML_d) e = cudaFree(tomo_gpu->thetaML_d);
  process_err(e, "free gpu thetaML_d");

  if (tomo_gpu->X_d) e = cudaFree(tomo_gpu->X_d);
  process_err(e, "free gpu X_d");

  if (tomo_gpu->Y_d) e = cudaFree(tomo_gpu->Y_d);
  process_err(e, "free gpu Y_d");

  if (tomo_gpu->XPup_d) e = cudaFree(tomo_gpu->XPup_d);
  process_err(e, "free gpu XPup_d");

  if (tomo_gpu->YPup_d) e = cudaFree(tomo_gpu->YPup_d);
  process_err(e, "free gpu YPup_d");

  /*
  if (tomo_gpu->Cmm_d) e = cudaFree(tomo_gpu->Cmm_d);
  process_err(e, "free gpu YPup_d");

  if (tomo_gpu->Cpm_d) e = cudaFree(tomo_gpu->Cpm_d);
  process_err(e, "free gpu YPup_d");

  if (tomo_gpu->R_d) e = cudaFree(tomo_gpu->R_d);
  process_err(e, "free gpu YPup_d");
  */

  if ((tomo_gpu->tabDPHI_d) != NULL) e = cudaFree(tomo_gpu->tabDPHI_d);
  process_err(e, "free gpu tabDPHI_d");

  if ((tomo_gpu->L0diff_d) != NULL) e = cudaFree(tomo_gpu->L0diff_d);
  process_err(e, "free gpu L0diff_d");

  // destroy matcov stream
  e = cudaStreamDestroy(tomo_gpu->matcov_stream);
  process_err(e, "destroy matcov stream");
}

//============================================================================================
//============================ CPHIM DPHI FUNCTIONS
//=============================
//============================================================================================
void tab_dphi_lowpass(double *tab_dphi, struct cphim_struct *cphim_struct,
                      long Ndphi, double *L0diff_d, int Nl0, double convert,
                      double convert_int)
// void tabulateDPHI_gpu_gb(double* tabDPHI_d, double* rr_d,struct tomo_struct
// tomo, long Ndphi, long *indexL0_h)
/* DOCUMENT tabDPHI = tabulateDPHI(rr,tomo,Ndphi, indexL0)
 <tomo>            :  structure with all the needed information
 <Ndphi>           :  size of rr
 <indexL0>         :  link between the index of the studied layer and the index
 of the precomputed one.

 Computes the phase structure function for a separation rr(x,y).
 The r0 is not taken into account : the final result of DPHI(x,y,L0)
 has to be scaled with r0^-5/3, with r0 expressed in meters, to get
 the right value.

 Computes the phase structure for each different L0 and give a array (indexL0)
 to link the index of the layer i and the index of tabDPHI : for the layer l,
 DPHI = DPHI( du, dv, indexL0[l],rr,tabDPHI, convert). SEE ALSO: DPHI
 */
{
  // Assume one thread per element
  int nblocks = (Ndphi * Nl0) / tabDPHI_thread_x +
                (((Ndphi * Nl0) % tabDPHI_thread_x) != 0);
  dim3 dimBlock(tabDPHI_thread_x, 1);
  dim3 dimGrid(nblocks, 1);

  // tabulateDPHI_lowpass_kernel<<<dimGrid, dimBlock, 0,
  // cphim_struct->cphim_stream>>>(tab_dphi,cphim_struct->tab_int_x,
  // cphim_struct->tab_int_y, cphim_struct->dx, L0diff_d, Nl0, Ndphi, convert,
  // convert_int, cphim_struct->int_npts);
  carmaCheckMsg("tabulateDPHI_gpu_gb_kernel<<<>>> execution failed\n");
  // CudaCheckError();
}

//============================================================================================
//============================ MATCOV V3/V4 DPHI/SUBAP FUNCTIONS
//=============================
//============================================================================================
void tab_dphi_gpu_gb(double *tab_dphi, struct gtomo_struct *tomo_gpu,
                     long Ndphi, double *L0diff_d, int Nl0, double convert)
// void tabulateDPHI_gpu_gb(double* tabDPHI_d, double* rr_d,struct tomo_struct
// tomo, long Ndphi, long *indexL0_h)
/* DOCUMENT tabDPHI = tabulateDPHI(rr,tomo,Ndphi, indexL0)
 <tomo>            :  structure with all the needed information
 <Ndphi>           :  size of rr
 <indexL0>         :  link between the index of the studied layer and the index
 of the precomputed one.

 Computes the phase structure function for a separation rr(x,y).
 The r0 is not taken into account : the final result of DPHI(x,y,L0)
 has to be scaled with r0^-5/3, with r0 expressed in meters, to get
 the right value.

 Computes the phase structure for each different L0 and give a array (indexL0)
 to link the index of the layer i and the index of tabDPHI : for the layer l,
 DPHI = DPHI( du, dv, indexL0[l],rr,tabDPHI, convert). SEE ALSO: DPHI
 */
{
  // Assume one thread per element
  int nblocks = (Ndphi * Nl0) / tabDPHI_thread_x +
                (((Ndphi * Nl0) % tabDPHI_thread_x) != 0);
  dim3 dimBlock(tabDPHI_thread_x, 1);
  dim3 dimGrid(nblocks, 1);

  tabulateDPHI_gpu_gb_kernel<<<dimGrid, dimBlock, 0, tomo_gpu->matcov_stream>>>(
      tab_dphi, L0diff_d, Nl0, Ndphi, convert);
  carmaCheckMsg("tabulateDPHI_gpu_gb_kernel<<<>>> execution failed\n");
  // CudaCheckError();
}
//------------------------------------------------------------------------------------
// extern "C"
void sub_pos_gpu_gb(struct gtomo_struct *tomo_gpu, long Nlayer)
// void subap_position_gpu_gb(struct tomo_struct tomo, double ***u, double ***v)
/* DOCUMENT DOCUMENT         subap_position(tomo, u, v)
   <tomo>                : structure with all the needed information.
   <u> and <v>           : 3d arrays containing the sub-apertures projected
   coordinates onto all the layers. u[0][2][1] is the X-coordinate of the subap
   2 of the WFS 0 on the layer 1.

   Computes the projected coordinates of all subapertures  projected onto all
   the layer
 */
{
  int msize = Nlayer * tomo_gpu->Nx;
  int nblocks = msize / tabDPHI_thread_x + ((msize % tabDPHI_thread_x) != 0);
  dim3 dimBlock(tabDPHI_thread_x, 1);
  dim3 dimGrid(nblocks, 1);
  /*
    int nb = (int)(2);
      long *tmp;
      tmp=(long*)malloc((nb)*sizeof(long));
      carmaSafeCall(cudaMemcpy(tmp, tomo_gpu->ioff_d, sizeof(long) * nb,
                    cudaMemcpyDeviceToHost));
      for (int ii = 0 ; ii < nb ; ii++){
          printf("%5.5d \n",tmp[ii]);
      }
  */
  // std::cout << "Nsubap : " << Nsubap << std::endl;
  subposition_gpu_gb_kernel<<<dimGrid, dimBlock, 0, tomo_gpu->matcov_stream>>>(
      tomo_gpu->Nw, tomo_gpu->Nx, tomo_gpu->Nsubap_d, Nlayer,
      tomo_gpu->alphaX_d, tomo_gpu->alphaY_d, tomo_gpu->h_d, tomo_gpu->GsAlt_d,
      tomo_gpu->Nssp_d, tomo_gpu->diamPup_d, tomo_gpu->thetaML_d,
      tomo_gpu->ioff_d, tomo_gpu->X_d, tomo_gpu->Y_d, tomo_gpu->XPup_d,
      tomo_gpu->YPup_d, tomo_gpu->u_d, tomo_gpu->v_d);
  carmaCheckMsg("subposition_gpu_gb_kernel<<<>>> execution failed\n");
  /*
     int nb = (int)tomo_gpu->Nx * Nlayer;
        double *tmpp;
        tmpp=(double*)malloc((nb)*sizeof(double));
        carmaSafeCall(cudaMemcpy(tmpp, tomo_gpu->v_d, sizeof(double) * nb,
                            cudaMemcpyDeviceToHost));
        for (int ii = 0 ; ii < nb ; ii++){
          printf("%5.5f \n",tmpp[ii]);
        }
        */
  // CudaCheckError();
}

//============================================================================================
//=============================== TOMO UPDATE FUNCTIONS
//======================================
//============================================================================================
void update_tomo_atm_gpu_gb(struct gtomo_struct *tomo_gpu,
                            sutra_sensors *sensors, sutra_atmos *atmos,
                            double *L0, double *cn2, double *alphaX,
                            double *alphaY) {
  cudaError_t e;

  double h[atmos->nscreens];
  int ii = 0;
  for (map<float, sutra_tscreen *>::iterator it = atmos->d_screens.begin();
       it != atmos->d_screens.end(); ++it) {
    h[ii] = (double)it->second->altitude;
    ii++;
  }
  // DEBUG_TRACE("Here !\n");
  double dmax = 0.0;
  double maxalt = h[atmos->nscreens - 1];
  long minssp = tomo_gpu->Nssp[0];
  for (int cc = 0; cc < tomo_gpu->Nw; cc++) {
    double tmp = sqrtf(alphaX[cc] * alphaX[cc] + alphaY[cc] * alphaY[cc]);
    if (tmp > dmax) dmax = tmp;
    if (minssp > tomo_gpu->Nssp[cc]) minssp = tomo_gpu->Nssp[cc];
  }
  const double crmax =
      dmax * 2 * maxalt + (1 + 1. / minssp) * tomo_gpu->DiamTel;
  const double pasDPHI = 1. / tomo_gpu->pasDPHI;  // inverse du pas de rr
  const long Ndphi = floor(crmax * pasDPHI) + 1;
  // const double convert = (double)(Ndphi-1)/(crmax+1./pasDPHI);

  e = cudaMemcpyAsync(tomo_gpu->h_d, h, atmos->nscreens * sizeof(double),
                      cudaMemcpyHostToDevice, tomo_gpu->matcov_stream);
  process_err(e, "copy gpu h_d");
  // DEBUG_TRACE("HERE !");

  e = cudaMemcpyAsync(tomo_gpu->cn2_d, cn2, atmos->nscreens * sizeof(double),
                      cudaMemcpyHostToDevice, tomo_gpu->matcov_stream);
  process_err(e, "copy gpu cn2_d");
  double *sspSizeL =
      (double *)malloc(sizeof(double) * tomo_gpu->Nw * atmos->nscreens);
  for (int cc = 0; cc < tomo_gpu->Nw * atmos->nscreens; cc++) {
    int n = cc / atmos->nscreens;
    int l = cc - n * atmos->nscreens;
    if (n >= sensors->nsensors()) n -= 1;
    sspSizeL[cc] = (((double)(tomo_gpu->DiamTel / sensors->d_wfs[n]->nxsub)) *
                    (1. - tomo_gpu->GsAlt[n] * h[l]));
  }
  // DEBUG_TRACE("HERE !");

  e = cudaMemcpyAsync(tomo_gpu->sspSizeL_d, sspSizeL,
                      tomo_gpu->Nw * atmos->nscreens * sizeof(double),
                      cudaMemcpyHostToDevice, tomo_gpu->matcov_stream);
  process_err(e, "copy gpu sspSizeL_d");
  cudaStreamSynchronize(tomo_gpu->matcov_stream);
  // Search the different L0 and build indexL0
  const long Nlayer = atmos->nscreens;
  long i, j;
  int cpt = 1;
  double tmp[Nlayer];
  long indexL0[Nlayer];
  tmp[0] = L0[0];
  indexL0[0] = 0;
  // DEBUG_TRACE("HERE !");

  for (i = 1; i < Nlayer; i++) {
    j = 0;
    const double l0 = L0[i];

    while ((j < cpt) && (tmp[j] != l0)) {
      j++;
    }

    indexL0[i] = j;

    if (j == cpt) {
      tmp[j] = l0;
      cpt++;
    }
  }
  e = cudaMemcpyAsync((tomo_gpu->indexL0_d), indexL0,
                      atmos->nscreens * sizeof(long), cudaMemcpyHostToDevice,
                      tomo_gpu->matcov_stream);
  process_err(e, "copy gpu indexL0_d");
  // DEBUG_TRACE("HERE !");

  int Nl0 = cpt;
  /*
  double L0diff[Nl0];
   //DEBUG_TRACE("Cpt = %d ",cpt);
  // allocate space for L0
  if ((tomo_gpu->L0diff_d) != NULL){cudaFree(tomo_gpu->L0diff_d);
  //DEBUG_TRACE("HERE !");

  e = cudaMalloc((void**)&(tomo_gpu->L0diff_d), Nlayer*sizeof(double));
  process_err(e, "alloc gpu L0diff_d");
  for (i = 0; i < Nl0; i++)  {
    L0diff[i] = tmp[i];
  }

  */
  if ((tomo_gpu->L0diff_d) != NULL) cudaFree(tomo_gpu->L0diff_d);
  e = cudaMalloc((void **)&(tomo_gpu->L0diff_d), Nlayer * sizeof(double));
  // DEBUG_TRACE("HERE !");

  // offload L0diff
  e = cudaMemcpyAsync(tomo_gpu->L0diff_d, L0, Nlayer * sizeof(double),
                      cudaMemcpyHostToDevice, tomo_gpu->matcov_stream);
  process_err(e, "offload L0diff");
  // précalcul de DPHI : que pour chaque différent L0
  if ((tomo_gpu->tabDPHI_d) != NULL) {
    cudaFree(tomo_gpu->tabDPHI_d);
  }
  // printf("tabDPHI alloc \n");
  e = cudaMalloc((void **)&(tomo_gpu->tabDPHI_d), Nl0 * Ndphi * sizeof(double));
  process_err(e, "alloc gpu tabDPHI_d");
  // DEBUG_TRACE("HERE !");
  // tab_dphi_gpu_gb(tomo_gpu->tabDPHI_d, tomo_gpu, Ndphi, tomo_gpu->L0diff_d,
  // Nl0,convert);
  // carmaSafeCall(cudaDeviceSynchronize());

  // %%%%%%% Computation of the sub-apertures positions and sizes %%%%%%%%%%%
  // u, v :arrays containing all the sub-apertures coordinates of all WFS, one
  // after the other u[0][1][3] is the X-coordinate of subap number 3 of wfs
  // number 0 at altitude 3
  // DEBUG_TRACE("HERE !");

  // Computes  u and v
  sub_pos_gpu_gb(tomo_gpu, (long)atmos->nscreens);
  // DEBUG_TRACE("HERE !");

  carmaSafeCall(cudaDeviceSynchronize());

  if (sspSizeL) free(sspSizeL);
  // DEBUG_TRACE("Here !\n");
}
//---------------------------------------------------------------------------------
void update_tomo_sys_gpu_gb(struct gtomo_struct *tomo_gpu,
                            sutra_sensors *sensors, double *alphaX,
                            double *alphaY) {
  cudaError_t e;

  long ioff[tomo_gpu->Nw];
  ioff[0] = 0;
  for (int i = 1; i < tomo_gpu->Nw; i++) {
    ioff[i] = ioff[i - 1] + sensors->d_wfs[i - 1]->nvalid;
  }

  e = cudaMemcpyAsync(tomo_gpu->ioff_d, ioff, tomo_gpu->Nw * sizeof(long),
                      cudaMemcpyHostToDevice, tomo_gpu->matcov_stream);
  process_err(e, "copy gpu ioff_d");

  e = cudaMemcpyAsync(tomo_gpu->alphaX_d, alphaX, tomo_gpu->Nw * sizeof(double),
                      cudaMemcpyHostToDevice, tomo_gpu->matcov_stream);
  process_err(e, "copy gpu alphaX_d");

  e = cudaMemcpyAsync(tomo_gpu->alphaY_d, alphaY, tomo_gpu->Nw * sizeof(double),
                      cudaMemcpyHostToDevice, tomo_gpu->matcov_stream);
  process_err(e, "copy gpu alphaY_d");

  e = cudaMemcpyAsync(tomo_gpu->GsAlt_d, tomo_gpu->GsAlt,
                      tomo_gpu->Nw * sizeof(double), cudaMemcpyHostToDevice,
                      tomo_gpu->matcov_stream);
  process_err(e, "copy gpu GsAlt_d");

  e = cudaMemcpyAsync(tomo_gpu->Nssp_d, tomo_gpu->Nssp,
                      tomo_gpu->Nw * sizeof(long), cudaMemcpyHostToDevice,
                      tomo_gpu->matcov_stream);
  process_err(e, "copy gpu Nssp_d");

  e = cudaMemcpyAsync(tomo_gpu->Nsubap_d, tomo_gpu->Nsubap,
                      tomo_gpu->Nw * sizeof(long), cudaMemcpyHostToDevice,
                      tomo_gpu->matcov_stream);
  process_err(e, "copy gpu Nsubap_d");

  e = cudaMemcpyAsync(tomo_gpu->diamPup_d, tomo_gpu->diamPup,
                      tomo_gpu->Nw * sizeof(double), cudaMemcpyHostToDevice,
                      tomo_gpu->matcov_stream);
  process_err(e, "copy gpu diamPup_d");

  e = cudaMemcpyAsync(tomo_gpu->XPup_d, tomo_gpu->XPup,
                      tomo_gpu->Nw * sizeof(double), cudaMemcpyHostToDevice,
                      tomo_gpu->matcov_stream);
  process_err(e, "copy gpu XPup_d");

  e = cudaMemcpyAsync(tomo_gpu->YPup_d, tomo_gpu->YPup,
                      tomo_gpu->Nw * sizeof(double), cudaMemcpyHostToDevice,
                      tomo_gpu->matcov_stream);
  process_err(e, "copy gpu YPup_d");

  e = cudaMemcpyAsync(tomo_gpu->thetaML_d, tomo_gpu->thetaML,
                      tomo_gpu->Nw * sizeof(double), cudaMemcpyHostToDevice,
                      tomo_gpu->matcov_stream);
  process_err(e, "copy gpu thetaML_d");
  // DEBUG_TRACE("Update \n");

  double *X;
  double *Y;
  int *tmpX;
  int *tmpY;
  X = (double *)malloc((tomo_gpu->Nx) * sizeof(double));
  Y = (double *)malloc((tomo_gpu->Nx) * sizeof(double));
  tmpX = (int *)malloc((tomo_gpu->Nx) * sizeof(int));
  tmpY = (int *)malloc((tomo_gpu->Nx) * sizeof(int));
  int ind = 0;
  double p2m;
  for (int i = 0; i < tomo_gpu->Nw; i++) {
    e = cudaMemcpyAsync(tmpX, sensors->d_wfs[i]->d_validsubsx->getData(),
                        sizeof(int) * sensors->d_wfs[i]->nvalid,
                        cudaMemcpyDeviceToHost, tomo_gpu->matcov_stream);
    process_err(e, "copy tmpX");
    e = cudaMemcpyAsync(tmpY, sensors->d_wfs[i]->d_validsubsy->getData(),
                        sizeof(int) * sensors->d_wfs[i]->nvalid,
                        cudaMemcpyDeviceToHost, tomo_gpu->matcov_stream);
    process_err(e, "copy tmpY");
    p2m = (tomo_gpu->DiamTel / (double)sensors->d_wfs[i]->nxsub) /
          ((double)(tmpX[1] - tmpX[0]));

    for (int j = 0; j < sensors->d_wfs[i]->nvalid; j++) {
      X[ind + j] = ((double)tmpX[j] * p2m) -
                   (double)((tomo_gpu->DiamTel / 2.) *
                            (1. - 1. / (double)sensors->d_wfs[i]->nxsub));
      Y[ind + j] = ((double)tmpY[j] * p2m) -
                   (double)((tomo_gpu->DiamTel / 2.) *
                            (1. - 1. / (double)sensors->d_wfs[i]->nxsub));
    }
    ind += sensors->d_wfs[i]->nvalid;
  }
  /*
    for (int ii = 0; ii<tomo_gpu->Nx ; ii++){
          std::cout << "X : " << X[ii] << std::endl;
    }
    for (int jj = 0; jj<tomo_gpu->Nx ; jj++){
          std::cout << "Y : " << Y[jj] << std::endl;
    }
  */
  // generateXY(tomo_gpu,sensors);
  e = cudaMemcpyAsync(tomo_gpu->X_d, X, tomo_gpu->Nx * sizeof(double),
                      cudaMemcpyHostToDevice, tomo_gpu->matcov_stream);
  process_err(e, "copy gpu X_d");
  e = cudaMemcpyAsync(tomo_gpu->Y_d, Y, tomo_gpu->Nx * sizeof(double),
                      cudaMemcpyHostToDevice, tomo_gpu->matcov_stream);
  process_err(e, "copy gpu Y_d");
  // cudaStreamSynchronize(tomo_gpu->matcov_stream);
  // DEBUG_TRACE("Update \n");
  cudaStreamSynchronize(tomo_gpu->matcov_stream);
  /*
    int nb = (int)(408);
    double *tmp;
    tmp=(double*)malloc((nb)*sizeof(double));
    carmaSafeCall(cudaMemcpy(tmp, tomo_gpu->Y_d, sizeof(double) * nb,
                    cudaMemcpyDeviceToHost));
    for (int ii = 0 ; ii < nb ; ii++){
          printf("%5.5f \n",tmp[ii]);
    }
  */
}

void generateXY(struct gtomo_struct *tomo, sutra_sensors *sensors)
/* DOCUMENT  generateXY(struct tomo_struct tomo, double *Nsubap)
 <tomo>               :  structure with all the needed information
 <tomo.X> & <tomo.Y>            :   arrays containing all the sub-apertures
coordinates of all WFS, one after the other <tomo.Nsubap>              :  number
of subaperture of ezach WFS Generate the position (X,Y) of each subapertures of
each WFS on the telescope pupil and the number of subaperture of ezach WFS
(Nsubap)
 */
{
  const double bornemin = -tomo->DiamTel / 2.;
  const double Rtel2 = (tomo->DiamTel * tomo->DiamTel) / 4.;
  long NsubapTot = 0;
  long n;

  // Total number of subapertures (without obstruction)
  for (n = 0; n < tomo->Nw; n++) {
    NsubapTot += tomo->Nssp[n] * tomo->Nssp[n];
  }
  // DEBUG_TRACE("Here !\n");
  const long cNsubapTot = NsubapTot;
  double x[cNsubapTot], y[cNsubapTot];
  int index[cNsubapTot];

  int cpt = 0;
  int ioff = 0;

  // Computation of all the subapertures' positions
  for (n = 0; n < tomo->Nw; n++) {
    long Nsap = 0;
    double pas = tomo->DiamTel / (1. * tomo->Nssp[n]);
    int i;
    double Robs2;

    // to avoid some bug that eliminates useful central subapertures when
    // obs=0.286
    if (tomo->Nssp[n] != 7 || (tomo->obs <= 0.285 || tomo->obs >= 0.29)) {
      Robs2 = tomo->DiamTel * tomo->obs / 2. * tomo->DiamTel * tomo->obs / 2.;
    } else {
      Robs2 = tomo->DiamTel * 0.285 / 2. * tomo->DiamTel * 0.285 / 2.;
    }
    // DEBUG_TRACE("Here !\n");
    if (tomo->Nssp[n] != 1) {
      for (i = 0; i < tomo->Nssp[n]; i++) {
        double tp =
            bornemin + pas / 2. * (2. * i + 1.);  // y-coord of current subap
        int j;

        for (j = 0; j < tomo->Nssp[n]; j++) {
          x[ioff + j] =
              bornemin + pas / 2. * (2. * j + 1.);  // x-coord of current subap
          y[ioff + j] = tp;

          double r2 = x[ioff + j] * x[ioff + j] + y[ioff + j] * y[ioff + j];
          // DEBUG_TRACE("Here !\n");
          // Search the non-valid subapertures
          if (r2 < Robs2 || r2 >= Rtel2) {
            index[cpt] = j + ioff;  // list of the useless subapertures index
            cpt++;
          } else {
            Nsap++;
          }
        }
        ioff += tomo->Nssp[n];
      }
      // tomo->Nsubap[n] = Nsap;
    } else {         // Special case (Nssp = 1)
      x[ioff] = 0.;  // x-coord of current subap
      y[ioff] = 0.;
      ioff += tomo->Nssp[n];
      // tomo->Nsubap[n] = 1;
    }
  }

  double *X;
  double *Y;
  std::cout << "sizeX = " << cNsubapTot - cpt << std::endl;
  X = (double *)malloc((cNsubapTot - cpt) * sizeof(double));
  Y = (double *)malloc((cNsubapTot - cpt) * sizeof(double));
  tomo->Nx = cNsubapTot - cpt;

  int a = 0;
  int off = 0;
  int borne = 0;
  int i;
  // Suppress the non-valid subapertures
  while (a <= cpt) {
    if (a == cpt) {
      borne = cNsubapTot;
    } else {
      borne = index[a];
    }

    for (i = off; i < borne; i++) {
      X[i - a] = x[i];
      Y[i - a] = y[i];
    }

    off = index[a] + 1;
    a++;
  }
  cudaError_t e;
  e = cudaMemcpyAsync(tomo->X_d, X, tomo->Nx * sizeof(double),
                      cudaMemcpyHostToDevice, tomo->matcov_stream);
  process_err(e, "copy gpu X_d");
  e = cudaMemcpyAsync(tomo->Y_d, Y, tomo->Nx * sizeof(double),
                      cudaMemcpyHostToDevice, tomo->matcov_stream);
  process_err(e, "copy gpu Y_d");
}

//============================================================================================
//============================= MATCOV 3
//=====================================================
//============================================================================================
// extern "C"
/*
void matcov_gpu_3(double* data, int nrows, int ncols, int xoffset, int yoffset,
int lda, struct tomo_struct tomo, struct gtomo_struct *tomo_gpu)
{
        /* *** matcov gpu kernel driver ***
        *  Arguments
        *  ==========
        *  data		double pointer: A pointer to the matrix/submatrix to be
generated. It
        *  			should always point to the first element in a
matrix/submatrix
        *
        *  nrows	integer: The number of rows of the matrix/submatrix to
be generated
        *
        *  ncols	integer: The number of columns of the matrix/submatrix
to be generated
        *
        *  xoffset	integer: The x-offset of the submatrix, must be zero if
the entire matrix
        *			is generated. Its the x-coordinate of the first
element in the matrix/submatrix
        *
        *  yoffset  integer: The y-offset of the submatrix, must be zero if the
entire matrix
        *			is generated. Its the y-coordinate of the first
element in the matrix/submatrix
        *
        *  lda		integer: The leading dimension of the matrix/submatrix
        */
/*
  //cudaError_t e;

  char uplo, copy;

  uplo = 'f';	// full generation is enabled by default
  copy = 'c';

  int type_mat = tomo.part;

  if(type_mat == 1) // Caa matrix
  {
        // check if a square diagonal tile is generated then we set uplo to 'l'
or 'u'
        // and then enable the copy
        // This also applies if the entire matrix will be generated
        // otherwise (off diagonal tile or non square submatrix) - full
generation is assumed
        if((xoffset == yoffset) && (nrows == ncols))	// if sqaure & diagonal
        {
                uplo = 'l';
                if(type_mat == 1)copy = 'c';
        }
        else	// full generation, copy is ignored
        {
                uplo = 'f';
        }
  }
  else if(type_mat == 2 || type_mat == 3) // Cmaa matrix
  {
        uplo = 'f';		// full generation, copy is ignored
  }
  else
  {
        printf("ERROR: unrecognized type_mat %d \n", type_mat); exit(1);
  }

  // %%%%%%% Pre-computation of DPHI %%%%%%%%%%
  //Computes an array of DPHI (tabDPHI) for an array of subaperture distance rr
for each DIFFERENT L0
  //const long Nw = tomo.Nw;
  const double crmax = tomo.rmax;
  const double pasDPHI = 1./tomo.pasDPHI; //inverse du pas de rr
  const long Ndphi = floor(crmax*pasDPHI)+1;
  const double convert = (double)(Ndphi-1)/(crmax+1./pasDPHI);

  //int size = tomo.Nslopes - 2 * tomo.Nsubap[tomo.Nw-1];

  int nbx = nrows / matcov_thread_x + (nrows%matcov_thread_x != 0);
  int nby = ncols / matcov_thread_y + (ncols%matcov_thread_y != 0);

  dim3 dimBlock(matcov_thread_x, matcov_thread_y);
  dim3 dimGrid(nbx, nby);
  const long Nsubap = tomo.Nsubap[0];

  // generate a full matrix
  matcov_kernel_3<<<dimGrid, dimBlock, 0, tomo_gpu->matcov_stream>>>(uplo, copy,
data, nrows, ncols, xoffset, yoffset, lda,
                                           convert,tomo_gpu->sspSizeL_d,tomo_gpu->Nssp_d,tomo_gpu->u_d,tomo_gpu->v_d,
                                           pasDPHI,tomo_gpu->tabDPHI_d,tomo_gpu->indexL0_d,tomo_gpu->cn2_d,
                                           Ndphi,tomo.Nw,atmos->nscreens,Nsubap,type_mat,tomo.DiamTel);

  //if (type_mat == 1)
  //  matcov_kernel_copy<<<dimGrid, dimBlock>>>(data, nrows, ncols, xoffset,
yoffset, lda);

  //cudaStreamSynchronize(tomo_gpu->matcov_stream);
}
*/

//============================================================================================
//=========================== MATCOV 4 (NOISE) KERNELS/FUNCTION
//==============================
//============================================================================================
__device__ double compute_element_4(
    int ipos, int jpos, double convert, double *sspSizeL, long *Nssp, double *u,
    double *v, double pasDPHI, double *tabDPHI, double *indexL0, double *cn2,
    int Ndphi, int Nw, int Nlayer, long *Nsubap_wfs, long Nx, double *alphaX,
    double *alphaY, double lgs_cst, double noise_var, double spotWidth,
    double dH_lgs, double alt_lgs, int type_mat, int nlgs, double teldiam) {
  /* *** Covariance matrix per-element generation ***
   *   Arguments
   *   =========
   *	ipos:		Integer: global x-coordinate of the element w.r.t. the
   *entire matrix jpos:		Integer: global y-coordinate of the element
   *w.r.t. the entire matrix
   */

  // for now return a dummy value

  const double lambda2 = 0.00026942094446267851;
  // long Nsubap = Nsubap_wfs[0];
  // WFS m

  long Nsubapx = Nsubap_wfs[0];
  int m = 0;
  while ((ipos / (2 * Nsubapx)) >= 1) {
    m++;
    Nsubapx += Nsubap_wfs[m];
  }
  Nsubapx -= Nsubap_wfs[m];

  // int m = ipos / (2 * Nsubap);
  if (type_mat == 3) m = Nw - 1;
  // WFS n

  long Nsubapy = Nsubap_wfs[0];
  int n = 0;
  while ((jpos / (2 * Nsubapy)) >= 1) {
    n++;
    Nsubapy += Nsubap_wfs[n];
  }
  Nsubapy -= Nsubap_wfs[n];

  // int n = jpos / (2 * Nsubap);
  if (type_mat == 2) n = Nw - 1;
  // subap i
  // int i = ipos % (2 * Nsubap_wfs[m]);
  int i = ipos - 2 * Nsubapx;
  // subap j
  // int j = jpos % (2 * Nsubap_wfs[n]);
  int j = jpos - 2 * Nsubapy;
  // xy i
  int xy_i;
  // xy j
  int xy_j;
  if (i >= Nsubap_wfs[m]) {
    i -= Nsubap_wfs[m];
    xy_i = 1;
  } else
    xy_i = 0;
  if (j >= Nsubap_wfs[n]) {
    j -= Nsubap_wfs[n];
    xy_j = 1;
  } else
    xy_j = 0;

  const double sspSizem = teldiam / Nssp[m];
  const double sspSizen = teldiam / Nssp[n];

  const double kk = lambda2 / (sspSizem * sspSizen);

  int type = xy_i * 2 + xy_j;

  // Layer l
  double covar = 0.0;
#pragma unroll
  for (int l = 0; l < Nlayer; l++) {
    double sspSizeml = sspSizeL[m * Nlayer + l];
    double sspSizenl = sspSizeL[n * Nlayer + l];
    // test if the altitude layers is not higher than the LGS altitude
    if ((sspSizeml > 0) && (sspSizenl > 0)) {
      int pos1 = i + Nsubapx + l * Nx;
      int pos2 = j + Nsubapy + l * Nx;
      // if(threadIdx.x == 6 && threadIdx.y == 0 && blockIdx.x == 6 &&
      // blockIdx.y == 1) if((pos1 >= 6840) || (pos2 >= 6839))
      //{
      //	printf("================ pos1 = %d, pos2 = %d \n", pos1, pos2);
      //}
      //(6,0,0) in block (0,2,0);
      double du = u[pos1] - u[pos2];
      double dv = v[pos1] - v[pos2];

      double s1 = sspSizeml * 0.5;
      double s2 = sspSizenl * 0.5;

      double ac = s1 - s2;
      double ad = s1 + s2;
      double bc = -ad;  // initially -s1-s2;
      double bd = -ac;  // initially -s1+s2;

      if (type == 0)
        covar += 0.5 /* pasDPHI*/ *
                 cov_XX_gpu_gb(du, dv, ac, ad, bc, bd, tabDPHI, indexL0[l],
                               convert, Ndphi) *
                 kk * cn2[l];
      else if (type == 3)
        covar += 0.5 /* pasDPHI*/ *
                 cov_YY_gpu_gb(du, dv, ac, ad, bc, bd, tabDPHI, indexL0[l],
                               convert, Ndphi) *
                 kk * cn2[l];
      else {  // if ((type == 1) || (type == 2))
        double s0 =
            sqrt(s1 * s1 + s2 * s2);  // half size of the subaperture equivalent
                                      // to a convolution by s1 and s2
        double dd =
            (s1 > s2) ? 1. - s2 / s1 : 1. - s1 / s2;  // Nono's style ....
        covar +=
            0.25 /* pasDPHI*/ *
            cov_XY_gpu_gb(du, dv, s0, tabDPHI, indexL0[l], convert, Ndphi) *
            kk * cn2[l] * (1. - dd * dd);
      }
    }
  }
  // adding noise

  if (m == n) {
    if (m < nlgs) {
      if (i == j) {
        // lgs case
        const int pos1 = i + Nsubapx;
        double x = u[pos1];
        double y = v[pos1];
        double xwfs = alphaX[m] * 206265;
        double ywfs = alphaY[m] * 206265;
        double lltx = 0;
        double llty = 0;
        const double lltnorm = sqrtf(xwfs * xwfs + ywfs * ywfs);
        if (lltnorm != 0) {
          lltx = xwfs / lltnorm * teldiam / 2.0;
          llty = ywfs / lltnorm * teldiam / 2.0;
        }
        x -= lltx;
        y -= llty;
        x = 206265. * dH_lgs * x / alt_lgs /
            alt_lgs;  // extension at Fwhm, in arcsec
        y = 206265. * dH_lgs * y / alt_lgs /
            alt_lgs;                           // extension at Fwhm, in arcsec
        double lgsExt = sqrtf(x * x + y * y);  // lengh of the extension
        double lgsTheta = x != 0 ? atanf(y / x) : 0.0;  // angle of extension
        double totalExt = sqrtf(lgsExt * lgsExt + spotWidth * spotWidth);
        // lengh of the extension including seeing, laser size, ...
        double ratio = totalExt / spotWidth;
        double noiseLongAxis = noise_var * ratio * ratio;
        if (type == 0)
          covar += noiseLongAxis * cosf(lgsTheta) * cosf(lgsTheta) +
                   noise_var * sinf(lgsTheta) * sinf(lgsTheta);
        else if (type == 3)
          covar += noiseLongAxis * sinf(lgsTheta) * sinf(lgsTheta) +
                   noise_var * cosf(lgsTheta) * cosf(lgsTheta);
        else
          covar +=
              (noiseLongAxis - noise_var) * sinf(lgsTheta) * cosf(lgsTheta);
      }
      if ((type == 0) || (type == 3)) covar += lgs_cst;
    } else {
      // ngs case
      if (i == j) {
        if ((type == 0) || (type == 3)) {
          covar += noise_var;
        }
      }
    }
  }

  return (double)covar;
}

__device__ double compute_cphim_element(
    int ipos, int jpos, double convert, double *sspSizeL, long *Nssp, double *u,
    double *v, double *xact, double *yact, double xref, double yref, long npts,
    double *L0, double *cn2, int Ndphi, int Nw, int Ndm, int Nlayer,
    long *Nsubap, long Nx, long *Nactu_tot, int Nact, long *NlayerDM,
    long *indLayerDm, double *dx, double *alphaX, double *alphaY,
    double lgs_cst, double noise_var, double spotWidth, double dH_lgs,
    double alt_lgs, double *Hlayer, double *Hdm, double FoV, int nlgs,
    double teldiam, double *k2, double *tab_int_x, double *tab_int_y) {
  /* *** Covariance matrix per-element generation ***
   *   Arguments
   *   =========
   *	ipos:		Integer: global x-coordinate of the element w.r.t. the
   *entire matrix jpos:		Integer: global y-coordinate of the element
   *w.r.t. the entire matrix
   */

  // for now return a dummy value

  const double lambda2 =
      0.016414031750058719;  // RASC * 0.5 * 1e-6 / 2. / pi ie lambda = 0.5e-6

  // DM m
  int m = 0;
  long Nactux = Nactu_tot[0];
  while ((ipos / Nactux) >= 1) {
    m++;
    Nactux += Nactu_tot[m];
  }
  Nactux -= Nactu_tot[m];
  // WFS n
  long Nsubapx = Nsubap[0];
  int n = 0;  // jpos / (2 * Nsubap);
  while ((jpos / (2 * Nsubapx)) >= 1) {
    n++;
    Nsubapx += Nsubap[n];
  }
  Nsubapx -= Nsubap[n];
  // if (type_mat == 2) n = Nw-1;

  // Nact i
  int i = ipos - Nactux;
  // subap j
  int j = jpos - 2 * Nsubapx;
  int type = 0;
  if (j >= Nsubap[n]) {
    j -= Nsubap[n];
    type = 1;
  }
  /*
  //xy i
  int xy_i;
  //xy j
  int xy_j;
  if (i>=Nact/Ndm) {
    i-= Nact/Ndm;
    xy_i = 1;
  } else xy_i = 0;
  if (j>=Nsubap) {
    j-= Nsubap;
    xy_j = 1;
  } else xy_j = 0;
  */
  const double sspSizen = teldiam / Nssp[n];

  const double kk =
      lambda2 * k2[m] / sspSizen;  // k2 = y_wfs(1).lambda / 2. / pi /
                                   // y_dm(1).unitpervolt (Yorick computation)

  // Layer l
  double covar = 0.0;
  // long Nlayer4Dm = NlayerDM[m];
  long otherDM = 0;
  for (int ll = 0; ll < m; ll++) otherDM += NlayerDM[ll];

#pragma unroll
  for (int l = 0; l < Nlayer; l++) {
    long dolayer = 0;

    if (indLayerDm[l] == m) dolayer = 1;

    if (dolayer) {
      double sspSizenl = sspSizeL[n * Nlayer + l];
      // test if the altitude layers is not higher than the LGS altitude
      if ((sspSizenl > 0)) {
        int pos_act = i + Nactux;  // + l * Nact;
        int pos_ssp = j + Nsubapx + l * Nx;
        // int pos_act = ipos;
        // int pos_ssp = jpos;
        double deltah = abs(Hlayer[l] - Hdm[m]);
        // double pDiam = teldiam + 2 * FoV * Hdm[m];
        // double hproj = pDiam / FoV;
        // double xproj = xact[pos_act]/hproj;
        // double yproj = yact[ppos_act]/hproj;
        // double dX = xproj*deltah;
        // double dY = yproj*deltah;
        double du = xact[pos_act] - u[pos_ssp];
        double dv = yact[pos_act] - v[pos_ssp];

        double s2 = sspSizenl * 0.5;
        double fc = 0.5 / sqrt(dx[m] * dx[m] + FoV * FoV * deltah * deltah);

        // double ac = 0.;
        // double ad = 2*s2;
        // double bc = -ad;   // initially -s1-s2;
        // double bd = 0;   // initially -s1+s2;

        if (type == 0)
          covar += 0.5 *
                   cphim_XX(du, dv, u[pos_ssp], v[pos_ssp], xref, yref, s2,
                            L0[l], fc, npts, tab_int_x, tab_int_y) *
                   kk * cn2[l];
        else
          covar += 0.5 *
                   cphim_YY(du, dv, u[pos_ssp], v[pos_ssp], xref, yref, s2,
                            L0[l], fc, npts, tab_int_x, tab_int_y) *
                   kk * cn2[l];
      }
    }
  }

  return covar;
}

//------------------------------------------------------------------------------------------
__global__ void matcov_kernel_4(char uplo, char copy, float *data, int nrows,
                                int ncols, int xoffset, int yoffset, int lda,
                                double convert, double *sspSizeL, long *Nssp,
                                double *u, double *v, double pasDPHI,
                                double *tabDPHI, double *indexL0, double *cn2,
                                int Ndphi, int Nw, int Nlayer, long *Nsubap,
                                long Nx, double *alphaX, double *alphaY,
                                double lgs_cst, double noise_var,
                                double spotWidth, double dH_lgs, double alt_lgs,
                                int type_mat, int nlgs, double teldiam) {
  /* *** covariance matrix generation kernel ***
   *	The kernel generates the element values in a given matrix/submatrix
   *   The generation function can be any function, as long as each element
   *   can be computed both individually and independently
   *
   *	see argument description in the kernel driver
   */

  // local thread coordinates w.r.t. thread block
  const int tx_ = threadIdx.x;
  const int ty_ = threadIdx.y;

  // local thread block coordinates w.r.t. kernel grid
  const int bx_ = blockIdx.x;
  const int by_ = blockIdx.y;

  // local coordinates of the element w.r.t. submatrix
  int lx = bx_ * blockDim.x + tx_;
  int ly = by_ * blockDim.y + ty_;

  // global coordinates of the elemnt w.r.t. the entire matrix
  int gx = lx + xoffset;
  int gy = ly + yoffset;

  // out-of-bound threads should terminate
  if ((lx >= nrows) || (ly >= ncols)) return;

  // Advance the data pointer accordingly
  // data += ly * lda + lx;

  double value;
  if (uplo == 'l') {
    if (gy <= gx) {
      value = compute_element_4(
          gx, gy, convert, sspSizeL, Nssp, u, v, pasDPHI, tabDPHI, indexL0, cn2,
          Ndphi, Nw, Nlayer, Nsubap, Nx, alphaX, alphaY, lgs_cst, noise_var,
          spotWidth, dH_lgs, alt_lgs, type_mat, nlgs, teldiam);
      data[ly * lda + lx] = (float)value;
      if (copy == 'c') data[lx * lda + ly] = (float)value;
    }
  } else if (uplo == 'u') {  // upper
    if (gx <= gy) {
      value = compute_element_4(
          gx, gy, convert, sspSizeL, Nssp, u, v, pasDPHI, tabDPHI, indexL0, cn2,
          Ndphi, Nw, Nlayer, Nsubap, Nx, alphaX, alphaY, lgs_cst, noise_var,
          spotWidth, dH_lgs, alt_lgs, type_mat, nlgs, teldiam);
      data[ly * lda + lx] = (float)value;
      if (copy == 'c') data[lx * lda + ly] = (float)value;
    }
  } else {  // uplo = 'f' full generation
    value = compute_element_4(gx, gy, convert, sspSizeL, Nssp, u, v, pasDPHI,
                              tabDPHI, indexL0, cn2, Ndphi, Nw, Nlayer, Nsubap,
                              Nx, alphaX, alphaY, lgs_cst, noise_var, spotWidth,
                              dH_lgs, alt_lgs, type_mat, nlgs, teldiam);
    data[ly * lda + lx] = (float)value;
  }

  // if ((type_mat == 3) || (gx <= gy)) {
  //  // call the generation function
  //  data[0] = compute_element_4(gx, gy, convert, sspSizeL, Nssp, u, v,
  //  pasDPHI, tabDPHI, indexL0, cn2, Ndphi, Nw, Nlayer,
  //				Nsubap, alphaX, alphaY, lgs_cst, noise_var,
  // spotWidth, dH_lgs, alt_lgs, type_mat, nlgs, teldiam); printf("gx = %d, gy =
  // %d -----
  // %.2f \n", gx, gy, data[0]);
  //}
}

__global__ void CPHIM_kernel(
    float *data, int nrows, int ncols, int xoffset, int yoffset, int lda,
    double convert, double *sspSizeL, long *Nssp, double *u, double *v,
    double *xact, double *yact, double xref, double yref, long npts, double *L0,
    double *cn2, int Ndphi, int Nw, int Ndm, int Nlayer, long *Nsubap, long Nx,
    long *Nactu, int Nact, long *NlayerDM, long *indLayerDm, double *dx,
    double *alphaX, double *alphaY, double lgs_cst, double noise_var,
    double spotWidth, double dH_lgs, double alt_lgs, double *Hlayer,
    double *Hdm, double FoV, int nlgs, double teldiam, double *k2,
    double *tab_int_x, double *tab_int_y) {
  /* *** covariance matrix generation kernel ***
   *	The kernel generates the element values in a given matrix/submatrix
   *   The generation function can be any function, as long as each element
   *   can be computed both individually and independently
   *
   *	see argument description in the kernel driver
   */

  // local thread coordinates w.r.t. thread block
  const int tx_ = threadIdx.x;
  const int ty_ = threadIdx.y;

  // local thread block coordinates w.r.t. kernel grid
  const int bx_ = blockIdx.x;
  const int by_ = blockIdx.y;

  // local coordinates of the element w.r.t. submatrix
  int lx = bx_ * blockDim.x + tx_;
  int ly = by_ * blockDim.y + ty_;

  // global coordinates of the elemnt w.r.t. the entire matrix
  int gx = lx + xoffset;
  int gy = ly + yoffset;

  // out-of-bound threads should terminate
  if ((lx >= nrows) || (ly >= ncols)) return;

  // Advance the data pointer accordingly
  // data += ly * lda + lx;

  double value;

  value = compute_cphim_element(
      gx, gy, convert, sspSizeL, Nssp, u, v, xact, yact, xref, yref, npts, L0,
      cn2, Ndphi, Nw, Ndm, Nlayer, Nsubap, Nx, Nactu, Nact, NlayerDM,
      indLayerDm, dx, alphaX, alphaY, lgs_cst, noise_var, spotWidth, dH_lgs,
      alt_lgs, Hlayer, Hdm, FoV, nlgs, teldiam, k2, tab_int_x, tab_int_y);
  data[ly * lda + lx] = (float)value;
}
//============================================================================================
//============================= MATCOV 4 (NOISE)
//=============================================
//============================================================================================
void matcov_gpu_4(float *data, int nrows, int ncols, int xoffset, int yoffset,
                  int lda, struct gtomo_struct *tomo_gpu, sutra_atmos *atmos,
                  sutra_sensors *sensors, double *alphaX, double *alphaY) {
  /* *** matcov gpu kernel driver ***
   *  Arguments
   *  ==========
   *  data		double pointer: A pointer to the matrix/submatrix to be
   *generated. It should always point to the first element in a matrix/submatrix
   *
   *  nrows	integer: The number of rows of the matrix/submatrix to be
   *generated
   *
   *  ncols	integer: The number of columns of the matrix/submatrix to be
   *generated
   *
   *  xoffset	integer: The x-offset of the submatrix, must be zero if the
   *entire matrix is generated. Its the x-coordinate of the first element in the
   *matrix/submatrix
   *
   *  yoffset  integer: The y-offset of the submatrix, must be zero if the
   *entire matrix is generated. Its the y-coordinate of the first element in the
   *matrix/submatrix
   *
   *  lda		integer: The leading dimension of the matrix/submatrix
   */

  // cudaError_t e;
  char uplo, copy;

  uplo = 'f';  // full generation is enabled by default
  copy = 'c';

  int type_mat = 1;

  if (type_mat == 1) {  // Caa matrix
    // check if a square diagonal tile is generated then we set uplo to 'l' or
    // 'u' and then enable the copy This also applies if the entire matrix will
    // be generated otherwise (off diagonal tile or non square submatrix) - full
    // generation is assumed
    if ((xoffset == yoffset) && (nrows == ncols)) {  // if sqaure & diagonal
      uplo = 'l';
      copy = 'c';
    } else {  // full generation, copy is ignored
      uplo = 'f';
    }
  }
  // else if(type_mat == 2) //
  else if (type_mat == 2 || type_mat == 3) {  // Cmaa matrix
    uplo = 'f';  // full generation, copy is ignored
  } else {
    printf("ERROR: unrecognized type_mat %d \n", type_mat);
    exit(1);
  }
  // uplo = 'f';
  // %%%%%%% Pre-computation of DPHI %%%%%%%%%%
  // Computes an array of DPHI (tabDPHI) for an array of subaperture distance rr
  // for each DIFFERENT L0
  double h[atmos->nscreens];
  int ii = 0;
  for (map<float, sutra_tscreen *>::iterator it = atmos->d_screens.begin();
       it != atmos->d_screens.end(); ++it) {
    h[ii] = (double)it->second->altitude;
    ii++;
  }

  double dmax = 0.0;
  double maxalt = h[atmos->nscreens - 1];
  int minssp = tomo_gpu->Nssp[0];
  for (int cc = 0; cc < tomo_gpu->Nw; cc++) {
    double tmp = sqrtf(alphaX[cc] * alphaX[cc] + alphaY[cc] * alphaY[cc]);
    if (tmp > dmax) dmax = tmp;
    if (tomo_gpu->Nssp[cc] < minssp) minssp = tomo_gpu->Nssp[cc];
  }
  const double crmax =
      dmax * 2 * maxalt + (1 + 1. / minssp) * tomo_gpu->DiamTel;

  const double pasDPHI = 1. / tomo_gpu->pasDPHI;  // inverse du pas de rr
  const long Ndphi = floor(crmax * pasDPHI) + 1;
  const double convert = (double)(Ndphi - 1) / (crmax + 1. / pasDPHI);

  int nbx = nrows / matcov_thread_x + (nrows % matcov_thread_x != 0);
  int nby = ncols / matcov_thread_y + (ncols % matcov_thread_y != 0);

  dim3 dimBlock(matcov_thread_x, matcov_thread_y);
  dim3 dimGrid(nbx, nby);
  const long Nsubap = sensors->d_wfs[0]->nvalid;  // tomo_gpu->Nx;
  /*
    int nb = (int)(atmos->nscreens*tomo_gpu->Nw);
    double *tmp;
    tmp=(double*)malloc((nb)*sizeof(double));
    carmaSafeCall(cudaMemcpy(tmp, tomo_gpu->sspSizeL_d, sizeof(double) * nb,
                    cudaMemcpyDeviceToHost));
    for (int ii = 0 ; ii < nb ; ii++){
          printf("%f \n",tmp[ii]);
    }
    //printf("convert : %9.9f\n",convert);
     */

  matcov_kernel_4<<<dimGrid, dimBlock, 0, tomo_gpu->matcov_stream>>>(
      uplo, copy, data, nrows, ncols, xoffset, yoffset, lda, convert,
      tomo_gpu->sspSizeL_d, tomo_gpu->Nssp_d, tomo_gpu->u_d, tomo_gpu->v_d,
      pasDPHI, tomo_gpu->tabDPHI_d, tomo_gpu->L0diff_d, tomo_gpu->cn2_d, Ndphi,
      tomo_gpu->Nw, atmos->nscreens, tomo_gpu->Nsubap_d, tomo_gpu->Nx,
      tomo_gpu->alphaX_d, tomo_gpu->alphaY_d, tomo_gpu->lgs_cst,
      (double)0.0 /*sensors->d_wfs[0]->noise*/, tomo_gpu->spot_width,
      tomo_gpu->lgs_depth, tomo_gpu->lgs_alt, type_mat, tomo_gpu->nlgs,
      tomo_gpu->DiamTel);
  carmaCheckMsg("matcov_kernel_4<<<>>> execution failed\n");
  cudaStreamSynchronize(tomo_gpu->matcov_stream);
  /*
  int nb = (int)sensors->d_wfs[0]->nvalid * 2;
  nb = nb*nb;
    double *tmp;
    tmp=(double*)malloc((nb)*sizeof(double));
    carmaSafeCall(cudaMemcpy(tmp, data, sizeof(double) * nb,
                    cudaMemcpyDeviceToHost));
    for (int ii = 0 ; ii < nb ; ii++)
        std::cout << tmp[ii] << std::endl;
        */
  // if (type_mat == 1)
  // matcov_kernel_copy<<<dimGrid, dimBlock>>>(data, nrows, ncols, xoffset,
  // yoffset, lda);
}

//============================================================================================
//============================= CPHIM
//=============================================
//============================================================================================
void CPHIM(float *data, int nrows, int ncols, int xoffset, int yoffset, int lda,
           struct cphim_struct *cphim_struct, sutra_atmos *atmos,
           sutra_sensors *sensors, double *alphaX, double *alphaY,
           carma_device *device) {
  /* *** matcov gpu kernel driver ***
   *  Arguments
   *  ==========
   *  data		double pointer: A pointer to the matrix/submatrix to be
   *generated. It should always point to the first element in a matrix/submatrix
   *
   *  nrows	integer: The number of rows of the matrix/submatrix to be
   *generated
   *
   *  ncols	integer: The number of columns of the matrix/submatrix to be
   *generated
   *
   *  xoffset	integer: The x-offset of the submatrix, must be zero if the
   *entire matrix is generated. Its the x-coordinate of the first element in the
   *matrix/submatrix
   *
   *  yoffset  integer: The y-offset of the submatrix, must be zero if the
   *entire matrix is generated. Its the y-coordinate of the first element in the
   *matrix/submatrix
   *
   *  lda		integer: The leading dimension of the matrix/submatrix
   */

  // %%%%%%% Pre-computation of DPHI %%%%%%%%%%
  // Computes an array of DPHI (tabDPHI) for an array of subaperture distance rr
  // for each DIFFERENT L0

  double h[atmos->nscreens];
  int ii = 0;
  for (map<float, sutra_tscreen *>::iterator it = atmos->d_screens.begin();
       it != atmos->d_screens.end(); ++it) {
    h[ii] = (double)it->second->altitude;
    ii++;
  }

  double dmax = 0.0;
  double maxalt = h[atmos->nscreens - 1];
  int minssp = cphim_struct->Nssp[0];
  for (int cc = 0; cc < cphim_struct->Nw; cc++) {
    double tmp = sqrtf(alphaX[cc] * alphaX[cc] + alphaY[cc] * alphaY[cc]);
    if (tmp > dmax) dmax = tmp;
    if (cphim_struct->Nssp[cc] < minssp) minssp = cphim_struct->Nssp[cc];
  }
  const double crmax =
      dmax * 2 * maxalt + (1 + 1. / minssp) * cphim_struct->DiamTel;

  const double pasDPHI = 1. / cphim_struct->pasDPHI;  // inverse du pas de rr
  const long Ndphi = floor(crmax * pasDPHI) + 1;
  const double convert = (double)(Ndphi - 1) / (crmax + 1. / pasDPHI);

  int nbx = nrows / matcov_thread_x + (nrows % matcov_thread_x != 0);
  int nby = ncols / matcov_thread_y + (ncols % matcov_thread_y != 0);

  dim3 dimBlock(matcov_thread_x, matcov_thread_y);
  dim3 dimGrid(nbx, nby);
  const long Nsubap = sensors->d_wfs[0]->nvalid;

  /*
    int nb = (int)(1224);
    double *tmp;
    tmp=(double*)malloc((nb)*sizeof(double));
    carmaSafeCall(cudaMemcpy(tmp, tomo_gpu->u_d, sizeof(double) * nb,
                    cudaMemcpyDeviceToHost));
    for (int ii = 0 ; ii < nb ; ii++){
          printf("%5.20f \n",tmp[ii]);
    }
    printf("convert : %9.9f\n",convert);
  */

  CPHIM_kernel<<<dimGrid, dimBlock, 0, cphim_struct->cphim_stream>>>(
      data, nrows, ncols, xoffset, yoffset, lda, convert,
      cphim_struct->sspSizeL_d, cphim_struct->Nssp_d, cphim_struct->u_d,
      cphim_struct->v_d, cphim_struct->xact_d, cphim_struct->yact_d,
      cphim_struct->x0, cphim_struct->y0, cphim_struct->int_npts,
      cphim_struct->L0diff_d, cphim_struct->cn2_d, Ndphi, cphim_struct->Nw,
      cphim_struct->Ndm, atmos->nscreens, cphim_struct->Nsubap_d,
      cphim_struct->Nx, cphim_struct->Nactu_tot_d, cphim_struct->Nactu,
      cphim_struct->NlayerDM_d, cphim_struct->indLayerDm_d, cphim_struct->dx_d,
      cphim_struct->alphaX_d, cphim_struct->alphaY_d, cphim_struct->lgs_cst,
      (double)0.0, cphim_struct->spot_width, cphim_struct->lgs_depth,
      cphim_struct->lgs_alt, cphim_struct->h_d, cphim_struct->hDm_d,
      cphim_struct->FoV, cphim_struct->nlgs, cphim_struct->DiamTel,
      cphim_struct->k2_d, cphim_struct->tab_int_x, cphim_struct->tab_int_y);

  carmaCheckMsg("matcov_kernel_4<<<>>> execution failed\n");
  cudaStreamSynchronize(cphim_struct->cphim_stream);
}

__device__ double unMoinsJ0(double x) {
  if (x < 0.1) {
    double x22 = (x * x) / 4.;
    return (1.0 - x22 / 4.) * x22;
  } else
    return (double)(1.0 - j0((double)x));
}

__global__ void compute_u831J0(double *x, double *y, int npts, double tmin,
                               double tmax, double dt) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  double t;
  while (tid < npts) {
    t = tmin + tid * dt;
    x[tid] = exp(t);
    y[tid] = exp(-t * (5 / 3.)) * unMoinsJ0(exp(t)) * dt;
    tid += blockDim.x * gridDim.x;
  }
}

__global__ void cuda_zcen_krnl(double *idata, double *odata, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    odata[tid] = (idata[tid + 1] + idata[tid]) / 2.0;
    tid += blockDim.x * gridDim.x;
  }
}

__global__ void intfrominftomin(double *data, double smallInt, int N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < N) {
    data[tid] += smallInt;
    tid += blockDim.x * gridDim.x;
  }
}

void sub_pos_cphim(struct cphim_struct *cphim_struct, long Nlayer)
// void subap_position_gpu_gb(struct tomo_struct tomo, double ***u, double ***v)
/* DOCUMENT DOCUMENT         subap_position(tomo, u, v)
   <tomo>                : structure with all the needed information.
   <u> and <v>           : 3d arrays containing the sub-apertures projected
   coordinates onto all the layers. u[0][2][1] is the X-coordinate of the subap
   2 of the WFS 0 on the layer 1.

   Computes the projected coordinates of all subapertures  projected onto all
   the layer
 */
{
  int msize = Nlayer * cphim_struct->Nx;
  int nblocks = msize / tabDPHI_thread_x + ((msize % tabDPHI_thread_x) != 0);
  dim3 dimBlock(tabDPHI_thread_x, 1);
  dim3 dimGrid(nblocks, 1);
  /*
    int nb = (int)(msize);
        double *tmpp;
        tmpp=(double*)malloc((nb)*sizeof(double));
        carmaSafeCall(cudaMemcpy(tmpp, cphim_struct->X_d, sizeof(double) * nb,
                            cudaMemcpyDeviceToHost));
        for (int ii = 0 ; ii < nb ; ii++){
          printf("%5.5f \n",tmpp[ii]);
        }
        */
  // std::cout << "Nsubap : " << Nsubap << std::endl;

  subposition_gpu_gb_kernel<<<dimGrid, dimBlock, 0,
                              cphim_struct->cphim_stream>>>(
      cphim_struct->Nw, cphim_struct->Nx, cphim_struct->Nsubap_d, Nlayer,
      cphim_struct->alphaX_d, cphim_struct->alphaY_d, cphim_struct->h_d,
      cphim_struct->GsAlt_d, cphim_struct->Nssp_d, cphim_struct->diamPup_d,
      cphim_struct->thetaML_d, cphim_struct->ioff_d, cphim_struct->X_d,
      cphim_struct->Y_d, cphim_struct->XPup_d, cphim_struct->YPup_d,
      cphim_struct->u_d, cphim_struct->v_d);

  carmaCheckMsg("subposition_gpu_gb_kernel<<<>>> execution failed\n");
  /*
     int nb = (int)(184);
        double *tmpp;
        tmpp=(double*)malloc((nb)*sizeof(double));
        carmaSafeCall(cudaMemcpy(tmpp, tomo_gpu->u_d, sizeof(double) * nb,
                            cudaMemcpyDeviceToHost));
        for (int ii = 0 ; ii < nb ; ii++){
          printf("%5.5f \n",tmpp[ii]);
        }
  */
  // CudaCheckError();
}

void init_cphim_struct(struct cphim_struct *cphim_struct, sutra_atmos *atmos,
                       sutra_sensors *sensors, sutra_dms *dms, double diamTel) {
  cudaError_t e;

  cphim_struct->DiamTel = diamTel;
  cphim_struct->Nw = sensors->nsensors();

  int Nactu = 0;
  int Ndm = 0;
  vector<sutra_dm *>::iterator p;
  p = dms->d_dms.begin();
  while (p != dms->d_dms.end()) {
    sutra_dm *dm = *p;
    if (dm->type != "tt") {
      Nactu += dm->ninflu;
      Ndm += 1;
    }
    p++;
  }

  cphim_struct->Nactu = Nactu;
  cphim_struct->Ndm = Ndm;
  cphim_struct->Nlayer = atmos->nscreens;
  cphim_struct->int_npts = 10000;
  cphim_struct->pasDu = 0.0001;
  cphim_struct->Nactu_tot = (long *)malloc(cphim_struct->Ndm * sizeof(long));
  cphim_struct->NlayerDM = (long *)malloc(cphim_struct->Ndm * sizeof(long));
  cphim_struct->indLayerDm = (long *)malloc(atmos->nscreens * sizeof(long));
  p = dms->d_dms.begin();
  int indx = 0;
  while (p != dms->d_dms.end()) {
    sutra_dm *dm = *p;
    if (dm->type != "tt") {
      cphim_struct->Nactu_tot[indx] = dm->ninflu;
      indx += 1;
    }
    p++;
  }

  cphim_struct->lgs_cst = 0.;
  cphim_struct->spot_width = 1.;
  cphim_struct->lgs_depth = 10000.;
  cphim_struct->lgs_alt = 90000.;
  cphim_struct->nlgs = 0;
  cphim_struct->pasDPHI = 0.0001;

  cphim_struct->Nx = 0;
  cphim_struct->Nssp = (long *)malloc(cphim_struct->Nw * sizeof(long));
  cphim_struct->Nsubap = (long *)malloc(cphim_struct->Nw * sizeof(long));
  cphim_struct->diamPup = (double *)malloc(cphim_struct->Nw * sizeof(double));
  cphim_struct->XPup = (double *)malloc(cphim_struct->Nw * sizeof(double));
  cphim_struct->YPup = (double *)malloc(cphim_struct->Nw * sizeof(double));
  cphim_struct->thetaML = (double *)malloc(cphim_struct->Nw * sizeof(double));
  cphim_struct->GsAlt = (double *)malloc(cphim_struct->Nw * sizeof(double));

  for (int i = 0; i < cphim_struct->Nw; i++) {
    cphim_struct->Nssp[i] = sensors->d_wfs[i]->nxsub;
    cphim_struct->Nsubap[i] = sensors->d_wfs[i]->nvalid;
    cphim_struct->diamPup[i] = (double)cphim_struct->Nssp[i];
    cphim_struct->XPup[i] = 0.;
    cphim_struct->YPup[i] = 0.;
    cphim_struct->thetaML[i] = 0.;
    if (sensors->d_wfs[i]->d_gs->lgs) {
      cphim_struct->nlgs += 1;
      cphim_struct->GsAlt[i] = 1.0 / cphim_struct->lgs_alt;
    } else
      cphim_struct->GsAlt[i] = 0.0;
    cphim_struct->GsAlt[i] = 0.;
    cphim_struct->Nx += sensors->d_wfs[i]->nvalid;
  }

  e = cudaMalloc((void **)&(cphim_struct->Nactu_tot_d),
                 cphim_struct->Ndm * sizeof(long));
  process_err(e, "alloc gpu Nactu_tot_d");

  e = cudaMalloc((void **)&(cphim_struct->NlayerDM_d),
                 cphim_struct->Ndm * sizeof(long));
  process_err(e, "alloc gpu NlayerDM_d");

  e = cudaMalloc((void **)&(cphim_struct->indLayerDm_d),
                 atmos->nscreens * sizeof(long));
  process_err(e, "alloc gpu Nactu_tot_d");

  e = cudaMalloc((void **)&(cphim_struct->indexL0_d),
                 atmos->nscreens * sizeof(long));
  process_err(e, "alloc gpu indexL0_d");

  e = cudaMalloc((void **)&(cphim_struct->u_d),
                 atmos->nscreens * cphim_struct->Nx * sizeof(double));
  process_err(e, "alloc gpu u_d");
  // printf("size of u is %d\n",atmos->nscreens * sensors->d_wfs[0]->nvalid *
  // cphim_struct->Nw); printf("u_d = 0x%x \n", (cphim_struct->u_d) );

  e = cudaMalloc((void **)&(cphim_struct->v_d),
                 atmos->nscreens * cphim_struct->Nx * sizeof(double));
  process_err(e, "alloc gpu v_d");
  // printf("size of v is %d\n", tomo.Nlayer*tomo.Nsubap[0]*tomo.Nw);
  // printf("v_d = 0x%x \n", (cphim_struct->v_d) );

  e = cudaMalloc((void **)&(cphim_struct->sspSizeL_d),
                 cphim_struct->Nw * atmos->nscreens * sizeof(double));
  process_err(e, "alloc gpu sspSizeL_d");

  e = cudaMalloc((void **)&(cphim_struct->cn2_d),
                 atmos->nscreens * sizeof(double));
  process_err(e, "alloc gpu cn2_d");

  e = cudaMalloc((void **)&(cphim_struct->h_d),
                 atmos->nscreens * sizeof(double));
  process_err(e, "alloc gpu h_d");

  e = cudaMalloc((void **)&(cphim_struct->hDm_d),
                 cphim_struct->Ndm * sizeof(double));
  process_err(e, "alloc gpu h_d");

  e = cudaMalloc((void **)&(cphim_struct->Nssp_d),
                 cphim_struct->Nw * sizeof(long));
  process_err(e, "alloc gpu Nssp_d");

  e = cudaMalloc((void **)&(cphim_struct->Nsubap_d),
                 cphim_struct->Nw * sizeof(long));
  process_err(e, "alloc gpu Nsubap_d");

  e = cudaMalloc((void **)&(cphim_struct->ioff_d),
                 cphim_struct->Nw * sizeof(long));
  process_err(e, "alloc gpu ioff_d");

  e = cudaMalloc((void **)&(cphim_struct->alphaX_d),
                 cphim_struct->Nw * sizeof(double));
  process_err(e, "alloc gpu alphaX_d");

  e = cudaMalloc((void **)&(cphim_struct->alphaY_d),
                 cphim_struct->Nw * sizeof(double));
  process_err(e, "alloc gpu alphaY_d");

  e = cudaMalloc((void **)&(cphim_struct->GsAlt_d),
                 cphim_struct->Nw * sizeof(double));
  process_err(e, "alloc gpu GsAlt_d");

  e = cudaMalloc((void **)&(cphim_struct->diamPup_d),
                 cphim_struct->Nw * sizeof(double));
  process_err(e, "alloc gpu diamPup_d");

  e = cudaMalloc((void **)&(cphim_struct->thetaML_d),
                 cphim_struct->Nw * sizeof(double));
  process_err(e, "alloc gpu thetaML_d");

  e = cudaMalloc((void **)&(cphim_struct->k2_d),
                 cphim_struct->Ndm * sizeof(double));
  process_err(e, "alloc gpu k2_d");

  e = cudaMalloc((void **)&(cphim_struct->X_d),
                 cphim_struct->Nx * sizeof(double));
  process_err(e, "alloc gpu X_d");
  // printf("size of X is %d\n", cphim_struct->Nx);
  // printf("X_d = 0x%x \n", (cphim_struct->X_d) );

  e = cudaMalloc((void **)&(cphim_struct->Y_d),
                 cphim_struct->Nx * sizeof(double));
  process_err(e, "alloc gpu Y_d");
  // printf("size of X is %d\n", tomo.Nx);
  // printf("Y_d = 0x%x \n", (cphim_struct->Y_d) );
  e = cudaMalloc((void **)&(cphim_struct->xact_d),
                 cphim_struct->Nactu * sizeof(double));
  process_err(e, "alloc gpu X_d");
  // printf("size of X is %d\n", cphim_struct->Nx);
  // printf("X_d = 0x%x \n", (cphim_struct->X_d) );

  e = cudaMalloc((void **)&(cphim_struct->yact_d),
                 cphim_struct->Nactu * sizeof(double));
  process_err(e, "alloc gpu Y_d");

  e = cudaMalloc((void **)&(cphim_struct->dx_d),
                 cphim_struct->Ndm * sizeof(double));
  process_err(e, "alloc gpu Y_d");

  e = cudaMalloc((void **)&(cphim_struct->XPup_d),
                 cphim_struct->Nw * sizeof(double));
  process_err(e, "alloc gpu XPup_d");

  e = cudaMalloc((void **)&(cphim_struct->YPup_d),
                 cphim_struct->Nw * sizeof(double));
  process_err(e, "alloc gpu YPup_d");

  cphim_struct->L0diff_d = NULL;
  cphim_struct->tabDPHI_d = NULL;

  e = cudaMalloc((void **)&(cphim_struct->tab_int_x),
                 cphim_struct->int_npts * sizeof(double));
  process_err(e, "alloc gpu tab_int_x");

  e = cudaMalloc((void **)&(cphim_struct->tab_int_y),
                 cphim_struct->int_npts * sizeof(double));
  process_err(e, "alloc gpu tab_int_y");

  e = cudaStreamCreate(&(cphim_struct->cphim_stream));
  process_err(e, "create cphim stream");
}

void tab_u831J0(double *tab_int_x, double *tab_int_y, long npts) {
  // DEBUG_TRACE("tab_int !\n");

  double tmin = -4.;
  double tmax = 10.;
  cudaError_t e;
  double *t = (double *)malloc(sizeof(double) * npts);
  e = cudaMemcpy(t, tab_int_x, (npts) * sizeof(double), cudaMemcpyDeviceToHost);
  process_err(e, "copy test");

  double *temp;
  temp = (double *)malloc((npts - 1) * sizeof(double));
  double *tab;
  tab = (double *)malloc((npts) * sizeof(double));
  double *temp_d;
  // int nblocks = 0 , nthreads = 0;
  double dt = (tmax - tmin) / (npts - 1);
  // getNumBlocksAndThreads(device, npts, nblocks, nthreads);
  /*
  int nblocks = msize / tabDPHI_thread_x + ( ( msize % tabDPHI_thread_x) != 0);
  dim3 dimBlock(tabDPHI_thread_x, 1);
  dim3 dimGrid(nblocks, 1);
  */
  int device;
  cudaGetDevice(&device);
  struct cudaDeviceProp props;
  cudaGetDeviceProperties(&props, device);
  int nthreads = props.maxThreadsPerBlock;
  int nblocks = npts / nthreads + ((npts % nthreads) != 0);
  dim3 grid(nblocks), threads(nthreads);
  compute_u831J0<<<grid, threads>>>(tab_int_x, tab_int_y, npts, tmin, tmax, dt);
  carmaCheckMsg("compute_u831J0<<<>>> execution failed\n");
  // DEBUG_TRACE("tab_int !\n");
  e = cudaMalloc((void **)&(temp_d), (npts - 1) * sizeof(double));
  process_err(e, "alloc gpu temp_d");
  nblocks = (npts - 1) / nthreads + (((npts - 1) % nthreads) != 0);
  dim3 grid2(nblocks);
  cuda_zcen_krnl<<<grid2, threads>>>(tab_int_y, temp_d, npts - 1);
  carmaCheckMsg("cuda_zcen_krnl<<<>>> execution failed\n");
  // cuda_zcen(tab_int_y,temp_d, npts-1, device);
  // DEBUG_TRACE("tab_int !\n");
  e = cudaMemcpy(temp, temp_d, (npts - 1) * sizeof(double),
                 cudaMemcpyDeviceToHost);
  process_err(e, "copy cpu temp");
  cudaFree(temp_d);
  cumsum(tab, temp, npts);
  // DEBUG_TRACE("tab_int !\n");
  e = cudaMemcpy(tab_int_y, tab, (npts) * sizeof(double),
                 cudaMemcpyHostToDevice);
  process_err(e, "copy gpu tab");
  // DEBUG_TRACE("tab_int !\n");
  double smallx = exp(tmin);
  double smallInt = (0.75 * pow(smallx, 1 / 3.) * (1 - smallx * smallx / 112.));
  // DEBUG_TRACE("tab_int !\n");
  intfrominftomin<<<grid, threads>>>(tab_int_y, smallInt, npts);
  carmaCheckMsg("intfrominftomin<<<>>> execution failed\n");
}

void cumsum(double *odata, double *idata, int N) {
  odata[0] = 0;
  for (int i = 1; i < N; i++) {
    odata[i] = idata[i - 1] + odata[i - 1];
  }
}

void cuda_zcen(double *idata, double *odata, int N, carma_device *device) {
  int nblocks = 0, nthreads = 0;
  getNumBlocksAndThreads(device, N, nblocks, nthreads);
  dim3 grid(nblocks), threads(nthreads);

  cuda_zcen_krnl<<<grid, threads>>>(idata, odata, N);
  carmaCheckMsg("cuda_zcen_krnl<<<>>> execution failed\n");
}

void update_cphim_atm(struct cphim_struct *cphim_struct, sutra_sensors *sensors,
                      sutra_atmos *atmos, double *L0, double *cn2,
                      double *alphaX, double *alphaY) {
  cudaError_t e;

  double h[atmos->nscreens];
  int ii = 0;
  for (map<float, sutra_tscreen *>::iterator it = atmos->d_screens.begin();
       it != atmos->d_screens.end(); ++it) {
    h[ii] = (double)it->second->altitude;
    ii++;
  }
  // DEBUG_TRACE("Here !\n");
  double dmax = 0.0;
  double maxalt = h[atmos->nscreens - 1];
  long minssp = cphim_struct->Nssp[0];
  for (int cc = 0; cc < cphim_struct->Nw; cc++) {
    double tmp = sqrtf(alphaX[cc] * alphaX[cc] + alphaY[cc] * alphaY[cc]);
    if (tmp > dmax) dmax = tmp;
    if (minssp > cphim_struct->Nssp[cc]) minssp = cphim_struct->Nssp[cc];
  }
  const double crmax =
      dmax * 2. * maxalt + (1 + 1. / minssp) * cphim_struct->DiamTel;
  const double pasDPHI = 1. / cphim_struct->pasDPHI;  // inverse du pas de rr
  const long Ndphi = floor(crmax * pasDPHI) + 1;
  cphim_struct->Ndphi = Ndphi;
  // const double convert = (double)(Ndphi-1)/(crmax+1./pasDPHI);
  // const double convert_int = (double)(cphim_struct->int_npts
  // -1)/(expf(10.0f)+expf(-(14./cphim_struct->int_npts))); const double
  // convert_int = 14./(cphim_struct->int_npts-1);
  // DEBUG_TRACE("Here !\n");
  e = cudaMemcpyAsync(cphim_struct->h_d, h, atmos->nscreens * sizeof(double),
                      cudaMemcpyHostToDevice, cphim_struct->cphim_stream);
  process_err(e, "copy gpu h_d");

  e = cudaMemcpyAsync(cphim_struct->cn2_d, cn2,
                      atmos->nscreens * sizeof(double), cudaMemcpyHostToDevice,
                      cphim_struct->cphim_stream);
  process_err(e, "copy gpu cn2_d");
  // DEBUG_TRACE("Here !\n");
  double *sspSizeL =
      (double *)malloc(sizeof(double) * cphim_struct->Nw * atmos->nscreens);
  for (int cc = 0; cc < cphim_struct->Nw * atmos->nscreens; cc++) {
    int n = cc / atmos->nscreens;
    int l = cc - n * atmos->nscreens;
    if (n >= sensors->nsensors()) n -= 1;
    sspSizeL[cc] =
        (((double)(cphim_struct->DiamTel / sensors->d_wfs[n]->nxsub)) *
         (1. -
          cphim_struct->GsAlt[n] *
              h[l]));  //+ 2*sqrt(alphaX[n]*alphaX[n]+alphaY[n]*alphaY[n])*h[l];
  }
  // DEBUG_TRACE("Here !\n");
  e = cudaMemcpyAsync(cphim_struct->sspSizeL_d, sspSizeL,
                      cphim_struct->Nw * atmos->nscreens * sizeof(double),
                      cudaMemcpyHostToDevice, cphim_struct->cphim_stream);
  process_err(e, "copy gpu sspSizeL_d");
  cudaStreamSynchronize(cphim_struct->cphim_stream);
  // Search the different L0 and build indexL0
  const long Nlayer = atmos->nscreens;
  long i, j;
  int cpt = 1;
  double tmp[Nlayer];
  long indexL0[Nlayer];
  tmp[0] = L0[0];
  indexL0[0] = 0;

  for (i = 1; i < Nlayer; i++) {
    j = 0;
    const double l0 = L0[i];

    while ((j < cpt) && (tmp[j] != l0)) {
      j++;
    }

    indexL0[i] = j;

    if (j == cpt) {
      tmp[j] = l0;
      cpt++;
    }
  }
  e = cudaMemcpyAsync((cphim_struct->indexL0_d), indexL0,
                      atmos->nscreens * sizeof(long), cudaMemcpyHostToDevice,
                      cphim_struct->cphim_stream);
  process_err(e, "copy gpu indexL0_d");
  int Nl0 = cpt;
  /*
  double L0diff[Nl0];
  // DEBUG_TRACE("Here !\n");
  // allocate space for L0
  process_err(e, "alloc gpu L0diff_d");
  for (i = 0; i < Nl0; i++)  {
    L0diff[i] = tmp[i];
  }
  */
  if ((cphim_struct->L0diff_d) != NULL) {
    cudaFree(cphim_struct->L0diff_d);
  }
  e = cudaMalloc((void **)&(cphim_struct->L0diff_d),
                 cphim_struct->Nlayer * sizeof(double));
  // offload L0diff
  e = cudaMemcpyAsync(cphim_struct->L0diff_d, L0,
                      cphim_struct->Nlayer * sizeof(double),
                      cudaMemcpyHostToDevice, cphim_struct->cphim_stream);
  process_err(e, "offload L0diff");
  // précalcul de DPHI : que pour chaque différent L0
  if ((cphim_struct->tabDPHI_d) != NULL) {
    cudaFree(cphim_struct->tabDPHI_d);
  }
  // printf("tabDPHI alloc \n");
  e = cudaMalloc((void **)&(cphim_struct->tabDPHI_d),
                 Nl0 * Ndphi * cphim_struct->Ndm * sizeof(double));
  process_err(e, "alloc gpu tabDPHI_d");

  // DEBUG_TRACE("%5.5d %5.5d %5.5f\n",Ndphi,Nl0,convert);
  /*
  int nb = (int)(3);
  // FILE *f = fopen("tabDPHI_d.txt","w");
    double *tmpp;
    tmpp=(double*)malloc((nb)*sizeof(double));
    carmaSafeCall(cudaMemcpy(tmpp, cphim_struct->h_d, sizeof(double) * nb,
                    cudaMemcpyDeviceToHost));
    for (int ii = 0 ; ii < nb ; ii++){
        DEBUG_TRACE("%5.5f \n",tmpp[ii]);
    }
    */

  // tab_dphi_lowpass(cphim_struct->tabDPHI_d, cphim_struct, Ndphi,
  // cphim_struct->L0diff_d, Nl0,convert,convert_int);
  // carmaSafeCall(cudaDeviceSynchronize());
  /*
  int nb = (int)(Ndphi);
   // FILE *f = fopen("tabDPHI_d.txt","w");
      double *tmpp;
      tmpp=(double*)malloc((nb)*sizeof(double));
      carmaSafeCall(cudaMemcpy(tmpp, cphim_struct->tabDPHI_d, sizeof(double) *
  nb, cudaMemcpyDeviceToHost)); for (int ii = 0 ; ii < nb ; ii++){ printf("%5.5f
  \n",tmpp[ii]);
      }
      */
  // %%%%%%% Computation of the sub-apertures positions and sizes %%%%%%%%%%%
  // u, v :arrays containing all the sub-apertures coordinates of all WFS, one
  // after the other u[0][1][3] is the X-coordinate of subap number 3 of wfs
  // number 0 at altitude 3

  // Computes  u and v
  // DEBUG_TRACE("Here %d %d %d!\n", (long)atmos->nscreens,
  // (long)cphim_struct->Nw, (long)sensors->d_wfs[0]->nvalid);
  sub_pos_cphim(cphim_struct, (long)atmos->nscreens);
  // carmaSafeCall(cudaDeviceSynchronize());
  /*
    int nb = (int)(atmos->nscreens * sensors->d_wfs[0]->nvalid *
    cphim_struct->Nw); double *tmpp; tmpp=(double*)malloc((nb)*sizeof(double));
      carmaSafeCall(cudaMemcpy(tmpp, cphim_struct->u_d, sizeof(double) * nb,
                    cudaMemcpyDeviceToHost));
      for (int ii = 0 ; ii < nb ; ii++){
          printf("%5.5f \n",tmpp[ii]);
      }
  */
  if (sspSizeL) free(sspSizeL);
  // DEBUG_TRACE("Here !\n");
}

void update_cphim_sys(struct cphim_struct *cphim_struct, sutra_sensors *sensors,
                      double *alphaX, double *alphaY, double *xactu,
                      double *yactu, double *X, double *Y, long *NlayerDm,
                      long *indLayerDm, double *alt_dm, double *pitch,
                      double *k2, double FoV) {
  cudaError_t e;

  long ioff[cphim_struct->Nw];
  ioff[0] = 0;
  for (int i = 1; i < cphim_struct->Nw; i++) {
    ioff[i] = ioff[i - 1] + sensors->d_wfs[i - 1]->nvalid;
  }

  cphim_struct->FoV = FoV;

  e = cudaMemcpyAsync(cphim_struct->Nactu_tot_d, cphim_struct->Nactu_tot,
                      cphim_struct->Ndm * sizeof(long), cudaMemcpyHostToDevice,
                      cphim_struct->cphim_stream);
  process_err(e, "copy gpu Nactu_tot_d");

  e = cudaMemcpyAsync(cphim_struct->hDm_d, alt_dm,
                      cphim_struct->Ndm * sizeof(double),
                      cudaMemcpyHostToDevice, cphim_struct->cphim_stream);
  process_err(e, "copy gpu hDm_d");

  e = cudaMemcpyAsync(cphim_struct->NlayerDM_d, NlayerDm,
                      cphim_struct->Ndm * sizeof(long), cudaMemcpyHostToDevice,
                      cphim_struct->cphim_stream);
  process_err(e, "copy gpu NlayerDM_d");

  e = cudaMemcpyAsync(cphim_struct->indLayerDm_d, indLayerDm,
                      cphim_struct->Nlayer * sizeof(long),
                      cudaMemcpyHostToDevice, cphim_struct->cphim_stream);
  process_err(e, "copy gpu indLayerDm_d");

  e = cudaMemcpyAsync(cphim_struct->k2_d, k2,
                      cphim_struct->Ndm * sizeof(double),
                      cudaMemcpyHostToDevice, cphim_struct->cphim_stream);
  process_err(e, "copy gpu k2_d");

  e = cudaMemcpyAsync(cphim_struct->ioff_d, ioff,
                      cphim_struct->Nw * sizeof(long), cudaMemcpyHostToDevice,
                      cphim_struct->cphim_stream);
  process_err(e, "copy gpu ioff_d");

  e = cudaMemcpyAsync(cphim_struct->alphaX_d, alphaX,
                      cphim_struct->Nw * sizeof(double), cudaMemcpyHostToDevice,
                      cphim_struct->cphim_stream);
  process_err(e, "copy gpu alphaX_d");

  e = cudaMemcpyAsync(cphim_struct->alphaY_d, alphaY,
                      cphim_struct->Nw * sizeof(double), cudaMemcpyHostToDevice,
                      cphim_struct->cphim_stream);
  process_err(e, "copy gpu alphaY_d");

  e = cudaMemcpyAsync(cphim_struct->GsAlt_d, cphim_struct->GsAlt,
                      cphim_struct->Nw * sizeof(double), cudaMemcpyHostToDevice,
                      cphim_struct->cphim_stream);
  process_err(e, "copy gpu GsAlt_d");

  e = cudaMemcpyAsync(cphim_struct->Nssp_d, cphim_struct->Nssp,
                      cphim_struct->Nw * sizeof(long), cudaMemcpyHostToDevice,
                      cphim_struct->cphim_stream);
  process_err(e, "copy gpu Nssp_d");

  e = cudaMemcpyAsync(cphim_struct->Nsubap_d, cphim_struct->Nsubap,
                      cphim_struct->Nw * sizeof(long), cudaMemcpyHostToDevice,
                      cphim_struct->cphim_stream);
  process_err(e, "copy gpu Nsubap_d");

  e = cudaMemcpyAsync(cphim_struct->diamPup_d, cphim_struct->diamPup,
                      cphim_struct->Nw * sizeof(double), cudaMemcpyHostToDevice,
                      cphim_struct->cphim_stream);
  process_err(e, "copy gpu diamPup_d");

  e = cudaMemcpyAsync(cphim_struct->XPup_d, cphim_struct->XPup,
                      cphim_struct->Nw * sizeof(double), cudaMemcpyHostToDevice,
                      cphim_struct->cphim_stream);
  process_err(e, "copy gpu XPup_d");

  e = cudaMemcpyAsync(cphim_struct->YPup_d, cphim_struct->YPup,
                      cphim_struct->Nw * sizeof(double), cudaMemcpyHostToDevice,
                      cphim_struct->cphim_stream);
  process_err(e, "copy gpu YPup_d");

  e = cudaMemcpyAsync(cphim_struct->thetaML_d, cphim_struct->thetaML,
                      cphim_struct->Nw * sizeof(double), cudaMemcpyHostToDevice,
                      cphim_struct->cphim_stream);
  process_err(e, "copy gpu thetaML_d");
  // DEBUG_TRACE("Update \n");
  /*
    double *X;
    double *Y;
    int *tmpX;
    int *tmpY;
    X=(double*)malloc((cphim_struct->Nx)*sizeof(double));
    Y=(double*)malloc((cphim_struct->Nx)*sizeof(double));
    tmpX=(int*)malloc((sensors->d_wfs[0]->nvalid)*sizeof(int));
    tmpY=(int*)malloc((sensors->d_wfs[0]->nvalid)*sizeof(int));
    int ind = 0;
    double p2m;
    for(int i=0 ; i<cphim_struct->Nw ; i++){
          if(i<sensors->nsensors()){
                  e =
    cudaMemcpyAsync(tmpX,sensors->d_wfs[i]->d_validsubsx->getData() ,
    sizeof(int) * sensors->d_wfs[i]->nvalid,
                    cudaMemcpyDeviceToHost,cphim_struct->cphim_stream);
                  process_err(e,"copy tmpX");
                  e =
    cudaMemcpyAsync(tmpY,sensors->d_wfs[i]->d_validsubsy->getData() ,
    sizeof(int) * sensors->d_wfs[i]->nvalid,
                            cudaMemcpyDeviceToHost,cphim_struct->cphim_stream);
                  process_err(e,"copy tmpY");
                  p2m =
    (cphim_struct->DiamTel/(double)sensors->d_wfs[i]->nxsub)/((double)(tmpX[1]-tmpX[0]));
          }
          else{
                  e =
    cudaMemcpyAsync(tmpX,sensors->d_wfs[i-1]->d_validsubsx->getData() ,
    sizeof(int) * sensors->d_wfs[i-1]->nvalid,
                                    cudaMemcpyDeviceToHost,cphim_struct->cphim_stream);
                  process_err(e,"copy tmpX");
                  e =
    cudaMemcpyAsync(tmpY,sensors->d_wfs[i-1]->d_validsubsy->getData() ,
    sizeof(int) * sensors->d_wfs[i-1]->nvalid,
                                            cudaMemcpyDeviceToHost,cphim_struct->cphim_stream);
                  process_err(e,"copy tmpY");
                  p2m =
    (cphim_struct->DiamTel/(double)sensors->d_wfs[i-1]->nxsub)/((double)(tmpX[1]-tmpX[0]));
          }

          for(int j=0 ; j<sensors->d_wfs[0]->nvalid ; j++){
                  if(i<sensors->nsensors()){
                          X[ind + j] = ((double)tmpX[j] * p2m)-
    (double)((cphim_struct->DiamTel/2.)*(1.-1./(double)sensors->d_wfs[i]->nxsub));
                          Y[ind + j] = ((double)tmpY[j] * p2m) -
    (double)((cphim_struct->DiamTel/2.)*(1.-1./(double)sensors->d_wfs[i]->nxsub));
                  }
                  else{
                          X[ind + j] = ((double)tmpX[j] * p2m)-
    (double)((cphim_struct->DiamTel/2.)*(1.-1./(double)sensors->d_wfs[i-1]->nxsub));
                          Y[ind + j] = ((double)tmpY[j] * p2m) -
    (double)((cphim_struct->DiamTel/2.)*(1.-1./(double)sensors->d_wfs[i-1]->nxsub));
                  }
          }
          if(i<sensors->nsensors())
          ind += sensors->d_wfs[i]->nvalid;
          else ind += sensors->d_wfs[i-1]->nvalid;
    }
    /*
    for (int ii = 0; ii<cphim_struct->Nx ; ii++){
          std::cout << "X : " << X[ii] << std::endl;
    }
    for (int jj = 0; jj<cphim_struct->Nx ; jj++){
          std::cout << "Y : " << Y[jj] << std::endl;
    }
    */
  // generateXY(cphim_struct,sensors);

  e = cudaMemcpyAsync(cphim_struct->X_d, X, cphim_struct->Nx * sizeof(double),
                      cudaMemcpyHostToDevice, cphim_struct->cphim_stream);
  process_err(e, "copy gpu X_d");
  e = cudaMemcpyAsync(cphim_struct->Y_d, Y, cphim_struct->Nx * sizeof(double),
                      cudaMemcpyHostToDevice, cphim_struct->cphim_stream);
  process_err(e, "copy gpu Y_d");

  e = cudaMemcpyAsync(cphim_struct->xact_d, xactu,
                      cphim_struct->Nactu * sizeof(double),
                      cudaMemcpyHostToDevice, cphim_struct->cphim_stream);
  process_err(e, "copy gpu xact_d");
  e = cudaMemcpyAsync(cphim_struct->yact_d, yactu,
                      cphim_struct->Nactu * sizeof(double),
                      cudaMemcpyHostToDevice, cphim_struct->cphim_stream);
  process_err(e, "copy gpu yact_d");
  e = cudaMemcpyAsync(cphim_struct->dx_d, pitch,
                      cphim_struct->Ndm * sizeof(double),
                      cudaMemcpyHostToDevice, cphim_struct->cphim_stream);
  process_err(e, "copy gpu dx_d");

  std::cout << "   computing tabulated integral...";
  tab_u831J0(cphim_struct->tab_int_x, cphim_struct->tab_int_y,
             cphim_struct->int_npts);
  std::cout << " done" << std::endl;
  cphim_struct->x0 = xactu[cphim_struct->Nactu / 2 + 1];
  cphim_struct->y0 = yactu[cphim_struct->Nactu / 2 + 1];
  // cphim_struct->dx = (xactu[1] - xactu[0]) * 0.5;
  // cudaStreamSynchronize(cphim_struct->cphim_stream);
  // DEBUG_TRACE("Update \n");
  cudaStreamSynchronize(cphim_struct->cphim_stream);
  /*
    int nb = (int)(cphim_struct->Nactu);
    double *tmp;
    tmp=(double*)malloc((nb)*sizeof(double));
    carmaSafeCall(cudaMemcpy(tmp, cphim_struct->yact_d, sizeof(double) * nb,
                    cudaMemcpyDeviceToHost));
    for (int ii = 0 ; ii < nb ; ii++){
          printf("%5.5f \n",tmp[ii]);
    }
  */
}

void free_cphim_struct(struct cphim_struct *cphim_struct) {
  cudaError_t e;

  if ((cphim_struct->dx_d)) e = cudaFree(cphim_struct->dx_d);
  process_err(e, "free gpu dx_d");

  if ((cphim_struct->hDm_d)) e = cudaFree(cphim_struct->hDm_d);
  process_err(e, "free gpu hDm_d");

  if ((cphim_struct->NlayerDM_d)) e = cudaFree(cphim_struct->NlayerDM_d);
  process_err(e, "free gpu NlayerDM_d");

  if ((cphim_struct->indLayerDm_d)) e = cudaFree(cphim_struct->indLayerDm_d);
  process_err(e, "free gpu indLayerDm_d");

  if ((cphim_struct->Nactu_tot_d)) e = cudaFree(cphim_struct->Nactu_tot_d);
  process_err(e, "free gpu Nactu_tot_d");

  if ((cphim_struct->u_d)) e = cudaFree(cphim_struct->u_d);
  process_err(e, "free gpu u_d");

  if (cphim_struct->v_d) e = cudaFree(cphim_struct->v_d);
  process_err(e, "free gpu v_d");

  if (cphim_struct->sspSizeL_d) e = cudaFree(cphim_struct->sspSizeL_d);
  process_err(e, "free gpu sspSizeL_d");

  if (cphim_struct->cn2_d) e = cudaFree(cphim_struct->cn2_d);
  process_err(e, "free gpu cn2_d");

  if (cphim_struct->h_d) e = cudaFree(cphim_struct->h_d);
  process_err(e, "free gpu h_d");

  if (cphim_struct->indexL0_d) e = cudaFree(cphim_struct->indexL0_d);
  process_err(e, "free gpu indexL0_d");

  if (cphim_struct->Nssp_d) e = cudaFree(cphim_struct->Nssp_d);
  process_err(e, "free gpu Nssp_d");

  if (cphim_struct->Nsubap_d) e = cudaFree(cphim_struct->Nsubap_d);
  process_err(e, "free gpu Nsubap_d");

  if (cphim_struct->ioff_d) e = cudaFree(cphim_struct->ioff_d);
  process_err(e, "free gpu ioff_d");

  if (cphim_struct->alphaX_d) e = cudaFree(cphim_struct->alphaX_d);
  process_err(e, "free gpu alphaX_d");

  if (cphim_struct->alphaY_d) e = cudaFree(cphim_struct->alphaY_d);
  process_err(e, "free gpu alphaY_d");

  if (cphim_struct->GsAlt_d) e = cudaFree(cphim_struct->GsAlt_d);
  process_err(e, "free gpu GsAlt_d");

  if (cphim_struct->diamPup_d) e = cudaFree(cphim_struct->diamPup_d);
  process_err(e, "free gpu diamPup_d");

  if (cphim_struct->thetaML_d) e = cudaFree(cphim_struct->thetaML_d);
  process_err(e, "free gpu thetaML_d");

  if (cphim_struct->X_d) e = cudaFree(cphim_struct->X_d);
  process_err(e, "free gpu X_d");

  if (cphim_struct->Y_d) e = cudaFree(cphim_struct->Y_d);
  process_err(e, "free gpu Y_d");

  if (cphim_struct->tab_int_x) e = cudaFree(cphim_struct->tab_int_x);
  process_err(e, "free gpu tab_int_x");

  if (cphim_struct->tab_int_y) e = cudaFree(cphim_struct->tab_int_y);
  process_err(e, "free gpu tab_int_y");

  if (cphim_struct->xact_d) e = cudaFree(cphim_struct->xact_d);
  process_err(e, "free gpu xact_d");

  if (cphim_struct->yact_d) e = cudaFree(cphim_struct->yact_d);
  process_err(e, "free gpu yact_d");

  if (cphim_struct->XPup_d) e = cudaFree(cphim_struct->XPup_d);
  process_err(e, "free gpu XPup_d");

  if (cphim_struct->YPup_d) e = cudaFree(cphim_struct->YPup_d);
  process_err(e, "free gpu YPup_d");

  /*
  if (cphim_struct->Cmm_d) e = cudaFree(cphim_struct->Cmm_d);
  process_err(e, "free gpu YPup_d");

  if (cphim_struct->Cpm_d) e = cudaFree(cphim_struct->Cpm_d);
  process_err(e, "free gpu YPup_d");

  if (cphim_struct->R_d) e = cudaFree(cphim_struct->R_d);
  process_err(e, "free gpu YPup_d");
  */

  if ((cphim_struct->tabDPHI_d) != NULL) e = cudaFree(cphim_struct->tabDPHI_d);
  process_err(e, "free gpu tabDPHI_d");

  if ((cphim_struct->L0diff_d) != NULL) e = cudaFree(cphim_struct->L0diff_d);
  process_err(e, "free gpu L0diff_d");

  // destroy matcov stream
  e = cudaStreamDestroy(cphim_struct->cphim_stream);
  process_err(e, "destroy matcov stream");
}

__global__ void test_dphi_highpass_krnl(double *odata, double *r, double *tabx,
                                        double *taby, double d, long N) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid < N) {
    double fc = 1 / (2 * d);
    odata[tid] = DPHI_highpass_gb(r[tid], fc, tabx, taby, N);
  }
}

void test_DPHI_highpass(double R, double x0, long npts, carma_device *device) {
  // DEBUG_TRACE("tab_int !\n");
  cudaError_t e;
  double *tabx_d;
  double *taby_d;
  e = cudaMalloc((void **)&(tabx_d), (npts) * sizeof(double));
  process_err(e, "copy tabx_d");
  e = cudaMalloc((void **)&(taby_d), (npts) * sizeof(double));
  process_err(e, "copy taby_d");
  tab_u831J0(tabx_d, taby_d, npts);

  double *tabx = (double *)malloc(npts * sizeof(double));
  double *taby = (double *)malloc(npts * sizeof(double));
  e = cudaMemcpy(tabx, tabx_d, (npts) * sizeof(double), cudaMemcpyDeviceToHost);
  process_err(e, "copy tabx");
  e = cudaMemcpy(taby, taby_d, (npts) * sizeof(double), cudaMemcpyDeviceToHost);
  process_err(e, "copy taby");
  FILE *myfile;
  myfile = fopen("tabulated_int.txt", "w");
  for (int i = 0; i < npts; i++) {
    fprintf(myfile, "%9.9f %9.9f\n", tabx[i], taby[i]);
  }
  fclose(myfile);

  double *odata;
  e = cudaMalloc((void **)&(odata), (npts) * sizeof(double));
  double *r;
  r = (double *)malloc(npts * sizeof(double));
  double dr = R / (npts - 1);
  for (int i = 0; i < npts; i++) {
    r[i] = i * dr;
  }
  double *r_d;
  e = cudaMalloc((void **)&(r_d), (npts) * sizeof(double));
  e = cudaMemcpy(r_d, r, (npts) * sizeof(double), cudaMemcpyHostToDevice);
  process_err(e, "copy cpu temp");
  process_err(e, "alloc gpu temp_d");
  int nblocks = 0;
  int nthreads = 0;
  getNumBlocksAndThreads(device, npts, nblocks, nthreads);
  dim3 grid(nblocks), threads(nthreads);
  test_dphi_highpass_krnl<<<grid, threads>>>(odata, r_d, tabx_d, taby_d, x0,
                                             npts);

  e = cudaMemcpy(tabx, r_d, (npts) * sizeof(double), cudaMemcpyDeviceToHost);
  process_err(e, "copy tabx");
  e = cudaMemcpy(taby, odata, (npts) * sizeof(double), cudaMemcpyDeviceToHost);
  process_err(e, "copy taby");
  myfile = fopen("DPHI_highpass.txt", "w");
  for (int i = 0; i < npts; i++) {
    fprintf(myfile, "%9.9f %9.9f\n", tabx[i], taby[i]);
  }
  fclose(myfile);

  cudaFree(tabx_d);
  cudaFree(taby_d);
  cudaFree(odata);
  cudaFree(r_d);
  free(r);
  free(tabx);
  free(taby);
}
