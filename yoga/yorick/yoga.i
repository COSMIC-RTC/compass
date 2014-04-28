//plug_dir, _("../","./",plug_dir());
plug_in,"yoga";
yoga_version = "1.0";

/*
  ___        _          ____ 
 / _ \ _ __ | |_   _   / ___|
| | | | '_ \| | | | | | |    
| |_| | | | | | |_| | | |___ 
 \___/|_| |_|_|\__, |  \____|
               |___/         

*/

//==================================================================

extern _dist
/* PROTOTYPE
   int _dist(pointer dptr, long dimx, long dimy, float xc, float yc)
*/

extern _poidev
/* PROTOTYPE
   void _poidev(float array xm, long n)
*/
/*
extern idx_test
/ * PROTOTYPE
   int idx_test(pointer a)
*/

extern snapTransformSize
/* PROTOTYPE
   int snapTransformSize(int dataSize)
 */

/*
  ____                           _   _   _ _   _ _ _ _   _           
 / ___| ___ _ __   ___ _ __ __ _| | | | | | |_(_) (_) |_(_) ___  ___ 
| |  _ / _ \ '_ \ / _ \ '__/ _` | | | | | | __| | | | __| |/ _ \/ __|
| |_| |  __/ | | |  __/ | | (_| | | | |_| | |_| | | | |_| |  __/\__ \
 \____|\___|_| |_|\___|_|  \__,_|_|  \___/ \__|_|_|_|\__|_|\___||___/
                                                                     
 */

//==================================================================
extern _yogaThreadExit
/* PROTOTYPE
  void _yogaThreadExit(void);
 */

extern _yogaThreadSync
/* PROTOTYPE
  void _yogaThreadSync(void);
 */

/*
quit_yorick = quit;
func quit(void) {
  _yogaThreadExit;
  quit_yorick;
}
*/

extern _yoga_init;
/* PROTOTYPE
   int _yoga_init(void)
*/
func yoga_init(void) {
  _yoga_init;
}

extern _yoga_start_profile;
/* PROTOTYPE
   int _yoga_start_profile(void)
*/
func yoga_start_profile(void) {
	_yoga_start_profile;
}

extern _yoga_stop_profile;
/* PROTOTYPE
   int _yoga_stop_profile(void)
*/
func yoga_stop_profile(void) {
	_yoga_stop_profile;
}

extern _yoga_setDevice;
/* PROTOTYPE
   int _yoga_setDevice(int)
*/
func yoga_setDevice(dev) {
  return _yoga_setDevice(dev);
}

extern _yoga_getDevice;
/* PROTOTYPE
   int _yoga_getDevice(void)
*/
func yoga_getDevice(void) {
  return _yoga_getDevice();
}

extern _yoga_getnDevice;
/* PROTOTYPE
   int _yoga_getnDevice(void)
*/
func yoga_getnDevice(void) {
  return _yoga_getnDevice();
}

extern activeDevice;
/* DOCUMENT activeDevice
   mydevice = activeDevice();
   activeDevice,mydevice;

   if called as a function returns the active cuda device
   if called as a subroutine, sets the active cuda device to the id
   given by the argument
     
   SEE ALSO:
 */

extern yoga_context;
/* DOCUMENT yoga_context
   context = yoga_context();

   init and stores the context
   
   SEE ALSO:
 */

/*
__   _____   ____    _       ___  _     _           _   
\ \ / / _ \ / ___|  / \     / _ \| |__ (_) ___  ___| |_ 
 \ V / | | | |  _  / _ \   | | | | '_ \| |/ _ \/ __| __|
  | || |_| | |_| |/ ___ \  | |_| | |_) | |  __/ (__| |_ 
  |_| \___/ \____/_/   \_\  \___/|_.__// |\___|\___|\__|
                                     |__/               
 */

extern yoga_obj;
/* DOCUMENT yoga_obj
   obj = yoga_obj("type",[ndims,dims1,dims2,...])
   or
   obj = yoga_obj(0.0f,[ndims,dims1,dims2,...])
   or
   obj = yoga_obj(my_array)
   or
   obj = yoga_obj(other_obj) 

   This function creates a new Yoga object either from a type and
   dimensions list, from a standard Yorick object or from another Yoga object
   1rst form: type is : float, double, etc ...
   2nd form: the first argument is a scalarof the same type as the object
   you want to create
   3rd form: my array is a standard Yorick array
   4rth form: other_obj is another yoga object
   
   SEE ALSO:
 */

extern yoga_host2device;
/* DOCUMENT yoga_host2device
   yoga_host2device,device_obj,host_data;

   This function transfers data in host_data to the pre-existing
   Yoga object device_obj
   
   SEE ALSO:
 */
extern yoga_device2host; 
/* DOCUMENT yoga_device2host
   host_data = yoga_device2host(device_obj [,opt_flag]);

   This function transfers data in Yoga object device_obj to Yorick
   array host_data. opt_flag (default 0) can be used to specify if the
   main data field (d_data, opt_flag unspecified or = 0) or the optional
   data field (o_data, opt_flag = 1) of the device object has to be
   transfered
        
   SEE ALSO:
 */

extern yoga_getp;
/* DOCUMENT yoga_getp
   ptr = yoga_getp(device_obj)

   This function creates a pointer of the device_obj
   
   SEE ALSO:
 */

extern yoga_getpl;
/* DOCUMENT yoga_getpl
   ptr = yoga_getpl(device_obj)

   This function creates a pointer of the device_obj returned as a long int
   
   SEE ALSO:
 */

extern yoga_setv;
/* DOCUMENT yoga_setv
   obj = yoga_setv(host_data)

   This function creates a cublasVector from the Yorick array host_data
   and stores it in new Yoga Object obj
   
   SEE ALSO:
 */
extern yoga_setm;
/* DOCUMENT yoga_setm
   obj = yoga_setm(host_data)

   This function creates a cublasMatrix from the Yorick array host_data
   and stores it in new Yoga Object obj
     
   SEE ALSO:
 */

/*
  extern yoga_getv;
  extern yoga_getm;
*/

/*cublas functions*/

extern yoga_imin;
/* DOCUMENT yoga_imin
   result = yoga_imin(obj)

   This function returns the smallest index of the element of the minimum
   magnitude (BLAS imin) of Yoga object obj
     
   SEE ALSO:
 */
extern yoga_imax;
/* DOCUMENT yoga_imax
   result = yoga_imax(obj)

   This function returns the smallest index of the element of the maximum
   magnitude (BLAS amax) of Yoga object obj
     
   SEE ALSO:
 */
extern yoga_asum;
/* DOCUMENT yoga_asum
   result = yoga_asum(obj)

   This function returns the sum of the absolute values (BLAS asum) of
   Yoga object obj
     
   SEE ALSO:
 */
extern yoga_sum;
/* DOCUMENT yoga_sum
   result = yoga_sum(obj)

   This function returns the sum of the Yoga object obj
     
   SEE ALSO:
 */
extern yoga_nrm2;
/* DOCUMENT yoga_nrm2
   result = yoga_nrm2(obj)

   This function returns the Euclidean norm (BLAS nrm2) of Yoga object obj
     
   SEE ALSO:
 */
extern yoga_scale;
/* DOCUMENT yoga_scale
   yoga_scale,obj,alpha

   This routine scales (BLAS scal) the vector in Yoga object obj by the
   scalar alpha. Can only be called as a subroutine.
     
   SEE ALSO:
 */
extern yoga_swap;
/* DOCUMENT yoga_swap
   yoga_swap,obj1,obj2

   This routine swaps (BLAS swap) the content of Yoga objects obj1 and obj2.
   Can only be called as a subroutine.
     
   SEE ALSO:
 */
extern yoga_copy;
/* DOCUMENT yoga_copy
   yoga_copy,obj1,obj2

   This routine copy the content of Yoga objects obj2 in obj1.
   Can only be called as a subroutine.
     
   SEE ALSO:
 */
extern yoga_axpy;
/* DOCUMENT yoga_axpy
   dest = yoga_axpy(src,alpha)
   or
   yoga_axpy,dest,alpha,src

   This routine multiplies vector src by alpha and adds it to vector dest
   (BLAS axpy). If it is called as a function, a new yoga object is created,
   if it is called as a subroutine, it works on pre-existing objects.
     
   SEE ALSO:
 */
extern yoga_axpy_cpu;
/* DOCUMENT yoga_axpyyoga_axpy_cpu
   dest = yoga_axpy(alpha, src1, src2)
   or
   yoga_axpy,dest,alpha,src

   This routine multiplies vector src by alpha and adds it to vector dest
   (BLAS axpy). If it is called as a function, a new yoga object is created,
   if it is called as a subroutine, it works on pre-existing objects.
     
   SEE ALSO:
 */
extern yoga_dot;
/* DOCUMENT yoga_dot
   res = yoga_dot(obj1,obj2)

   This function returns the dot product (BLAS dot) of Yoga objects obj1
   and obj2.
     
   SEE ALSO:
 */
extern yoga_mv;
/* DOCUMENT yoga_mv
   vecty = yoga_mv(matA,vectx[,alpha])
   or
   yoga_mv,vecty,matA,vectx[,alpha,beta]

   This function performs the matrix-vector multiplication product (BLAS gemv)
   y = alpha * A * x + beta * y
   If called as a function, it creates a new yoga object
     
   SEE ALSO:
 */
extern yoga_mv_sparse;
/* DOCUMENT yoga_mv_sparse
   vecty = yoga_mv_sparse(matA,vectx[,alpha])
   or
   yoga_mv_sparse,vecty,matA,vectx[,alpha,beta]

   This function performs the matrix-vector multiplication product (BLAS gemv)
   y = alpha * A * x + beta * y
   If called as a function, it creates a new yoga object
     
   SEE ALSO:
 */
extern yoga_symv;
/* DOCUMENT yoga_symv
   vecty = yoga_symv(matA,vectx[,alpha])
   or
   yoga_symv,vecty,matA,vectx[,alpha,beta]

   This function performs the matrix-vector multiplication product (BLAS symv)
   the matrix A is considered symmetric
   y = alpha * A * x + beta * y
   If called as a function, it creates a new yoga object
     
   SEE ALSO:
 */
extern yoga_rank1;
/* DOCUMENT yoga_rank1
   matA = yoga_rank1(vectx,vecty)
   or
   yoga_rank1,matA,vectx,vecty

   This function performs the rank-1 update (BLAS ger)
   A = x * tranpose(y) + A
   If called as a function, it creates a new yoga object
     
   SEE ALSO:
 */
extern yoga_mm;
/* DOCUMENT yoga_mm
   matC = yoga_mm(matA,matB [,opA] [,opB] [,alpha])
   or
   yoga_mm,matC,matA,matB [,opA] [,opB] [,alpha] [,beta] 

   This function performs the matrix-matrix multiplication (BLAS gemm)
   C = alpha * opA( A ) * opB ( B ) + beta * C

   opA : operation on matrix A : 't' for transpose, 'n' for nothing
   opB : operation on matrix B : 't' for transpose, 'n' for nothing

   if you need to perform an operation on B you need to provide the operation on A as well
   if you need to provide alpha, you need to provide operations on A and B (even if default 'n')
   if you need to provide beta, you need to provide all the previous arguments (opA, opB, alpha)
   
   If called as a function, it creates a new yoga object so the beta argument is not applicable
   
   SEE ALSO:
 */

extern yoga_mm_cpu;
/* DOCUMENT yoga_mm_cpu
   matC = yoga_mm_cpu(matA,matB [,opA] [,opB] [,alpha])
   or
   yoga_mm_cpu,matC,matA,matB [,opA] [,opB] [,alpha] [,beta] 

   This function performs the matrix-matrix multiplication (BLAS gemm)
   C = alpha * opA( A ) * opB ( B ) + beta * C

   opA : operation on matrix A : 't' for transpose, 'n' for nothing
   opB : operation on matrix B : 't' for transpose, 'n' for nothing

   if you need to perform an operation on B you need to provide the operation on A as well
   if you need to provide alpha, you need to provide operations on A and B (even if default 'n')
   if you need to provide beta, you need to provide all the previous arguments (opA, opB, alpha)
   
   If called as a function, it creates a new yoga object so the beta argument is not applicable
   
   SEE ALSO:
 */

extern yoga_symm;
/* DOCUMENT yoga_symm
   matC = yoga_symm(matA,matB [,side] [,alpha])
   or
   yoga_symm,matC,matA,matB [,side] [,alpha] [,beta] 

   This function performs the symmetric matrix-matrix multiplication (BLAS symm)
   C = alpha * A * B  + beta * C if side == [] or 'l'
   C = alpha * B * A  + beta * C if side == 'r'
   with A a symmetric matrix

   if you need to provide alpha, you need to provide side (even if default 'l')
   if you need to provide beta, you need to provide all the previous arguments (side, alpha)
   
   If called as a function, it creates a new yoga object so the beta argument is not applicable
   
   SEE ALSO:
 */

extern yoga_syrk;
/* DOCUMENT yoga_syrk
   matC = yoga_syrk(matA,[,opA] [,alpha])
   or
   yoga_syrk,matC,matA,[,opA] [,alpha] [,beta] 

   This function performs the symmetric matrix rank1 operation (BLAS syrk)
   C = alpha * opA( A ) * transpose( opA ( A )) + beta * C

   opA : operation on matrix A : 't' for transpose, 'n' for nothing

   if you need to provide alpha, you need to provide operation on A (even if default 'n')
   if you need to provide beta, you need to provide all the previous arguments (opA, alpha)
   
   If called as a function, it creates a new yoga object so the beta argument is not applicable
   
   SEE ALSO:
 */

extern yoga_syrkx;
/* DOCUMENT yoga_syrkx
   matC = yoga_syrkx(matA,matB [,op] [,alpha])
   or
   yoga_syrkx,matC,matA,matB [,opA] [,alpha] [,beta] 

   This function performs the matrix-matrix multiplication (BLAS syrkx)
   the resulting matrix C needs to be symmetric which impose conditions on B
   
   C = alpha * op( A ) * transpose( op ( B ) )+ beta * C

   op : operation on matrix : 't' for transpose, 'n' for nothing

   if you need to perform an operation on B you need to provide the operation on A as well
   if you need to provide alpha, you need to provide operation on A and B (even if default 'n')
   if you need to provide beta, you need to provide all the previous arguments (op, alpha)
   
   If called as a function, it creates a new yoga object so the beta argument is not applicable
   
   SEE ALSO:
 */

extern yoga_am;
/* DOCUMENT yoga_am
   matC = yoga_am(matA,matB [,opA] [,opB] [,alpha] [,beta])
   or
   yoga_am,matC,matA,matB [,opA] [,opB] [,alpha] [,beta] 

   This function performs the matrix-matrix addition 
   
   C = alpha * opA( A ) * beta * opB( B )

   op : operation on matrix : 't' for transpose, 'n' for nothing

   if you need to perform an operation on B you need to provide the operation on A as well
   if you need to provide alpha, you need to provide operation on A and B (even if default 'n')
   if you need to provide beta, you need to provide all the previous arguments (op, alpha)
   
   If called as a function, it creates a new yoga object
   
   SEE ALSO:
 */

extern yoga_dmm;
/* DOCUMENT yoga_dmm
   matC = yoga_am(matA,vectx [,side])
   or
   yoga_am,matC,matA,vectx [,side] 

   This function performs the matrix-diagonal matrix multiply
   diagonal matrix is provided through the diagonal vector X
   
   C =  A * diag(X) if side == 'l'
   C =  diag(X) * A if side == 'r'

   default is side = 'l'

   If called as a function, it creates a new yoga object
   
   SEE ALSO:
 */

/*custom transpose*/
extern yoga_transpose;
/* DOCUMENT yoga_transpose
   matB = yoga_transpose(matA)
   or
   yoga_transpose,matB, matA

   This function performs the tranpose of yoga object matA and stores it
   in matB
   If called as a function, it creates a new yoga object
     
   SEE ALSO:
 */

/*random numbers*/
extern yoga_random;
/* DOCUMENT yoga_random
   obj = yoga_random("type",[ndims,dims1,dims2,...])
   or
   obj = yoga_random(0.0f,[ndims,dims1,dims2,...])
   or
   yoga_random,obj

   This function fills yoga object obj with random numbers following a
   uniform distribution.
   If called as a function, it creates a new yoga object
     
   SEE ALSO:
 */
extern yoga_random_n;
/* DOCUMENT yoga_random_n
   obj = yoga_random_n("type",[ndims,dims1,dims2,...])
   or
   obj = yoga_random_n(0.0f,[ndims,dims1,dims2,...])
   or
   yoga_random_n,obj

   This function fills yoga object obj with random numbers following a
   normal distribution.
   If called as a function, it creates a new yoga object
     
   SEE ALSO:
 */

/*fourier transform*/
extern yoga_fft;
/* DOCUMENT yoga_fft
   dest = yoga_fft(src[,dir])
   or
   yoga_fft,src,dest[,dir]
   or
   yoga_fft,src

   this routine performs the fft of yoga object src
   if called as a subroutine it either performs an inplace transform or an out-of-place
   transform if dest is specified
   if called as a function, it performs an out-of-place transform and creates a new
   yoga object to store the result
   if src is a 3-dimensional array, this routines assumes that you want to perfrom
   a batched fft of multiple 2d arrays, the number being specified by the 3rd dim
   
   SEE ALSO: yoga_fft.i
 */


/*general utilities on arrays*/
extern yoga_getarray;
/* DOCUMENT yoga_getarray
   dest = yoga_getarray(src[,xmin:xmax,ymin:ymax])
   or
   yoga_getarray,dest,src[,xmin:xmax,ymin:ymax]

   This routine extracts a sub array in src as specified by optional
   ranges : xmin:xmax and ymin:ymax
   
   if called as a function it returns the subarray in a new yoga object
   if called as a subroutine it uses pre-existing arrays

   SEE ALSO: 
 */
extern yoga_fillarray;
/* DOCUMENT yoga_fillarray
   yoga_fillarray,dest,src[,xmin:xmax,ymin:ymax]

   This routine fills a sub array of dest as specified by optional
   ranges : xmin:xmax and ymin:ymax by the content of src

   can only be called as a subroutine
   

   SEE ALSO: 
 */
extern yoga_getvalue;
/* DOCUMENT yoga_getvalue
   res = yoga_getvalue(obj,position)

   This routine returns the value of obj at the index position

   can only be called as a function. useful for debugs

   SEE ALSO: 
 */
extern yoga_plus;
/* DOCUMENT yoga_plus
   yoga_plus,obj,alpha

   This routine adds scalar alpha to all the elements of obj

   can only be called as a subroutine

   SEE ALSO: 
 */
extern yoga_plusai;
/* DOCUMENT yoga_plusai
   yoga_plusai,dest,src,ind

   This routine adds scalar src[ind] to all the elements of dest

   can only be called as a subroutine

   SEE ALSO: 
 */

/*fft convolution*/
extern yoga_fftconv_init;
/* DOCUMENT yoga_fftconv_init
   pdata_obj = yoga_fftconv_init(im_obj,im_ker,"real")
   or
   pspec_obj = yoga_fftconv_init(im_obj,im_ker,"complex")

   This routine inits the workspace for fft convolution of im_obj by ker_obj
   im_obj and ker_obj can be of different size
   it returns either the object that will contain the padded data
   or the array that will contain the padded spectrum (both are needed)
   can only be called as a function

   SEE ALSO: yoga_fft.i
 */
extern yoga_fftconv;
/* DOCUMENT yoga_fftconv
   yoga_fftconv,res_obj,im_obj,im_ker,pdata_obj,pspec_obj,kerX,kerY

   This routine computes the fft convolution of im_obj by ker_obj
   and stores it in res_obj
   im_obj and ker_obj can be of different size
   requires pdata_obj and pspec_obj to be initialized using yoga_fftconv_init
   kerX and kerY give the position of the kernel center in the im_obj space.
   can only be called as a subroutine

   SEE ALSO: yoga_fft.i
 */

extern yoga_test;

extern yoga_svd;
/* DOCUMENT yoga_svd
   yoga_svd,mat,eigenvals,vt,u

   This function computes the svd of matrix mat
   mat   : m x n
   VT    : n x n
   U     : m x m
   eigen : min(m,n)
   
   SEE ALSO:
 */

extern yoga_syevd;
/* DOCUMENT yoga_syevd
      yoga_syevd,mat,eigenvals, U, noComputeU= 
   OR yoga_syevd,mat,eigenvals, noComputeU=
   
   This function computes the svd of matrix mat
   mat   : n x n (CArMA obj) (WARNING:change into U if U not given)
   eigen : 1 x n (yArray)
   U     : n x n (CArMA obj) (optional)

   Even if noComputeU flag is set to 1, the matrix mat will be changed if U is not specified.

   SEE ALSO:
 */

extern yoga_syevd_m;
/* DOCUMENT yoga_syevd_m
   yoga_syevd_m, ngpu,mat,eigenvals,u, noComputeU= 

   This function computes the svd of matrix mat
   ngpu  : number of gpu
   mat   : n x n (yArray)
   eigen : 1 x n (yArray)
   U     : n x n (yArray)
   
   SEE ALSO:
 */

extern yoga_getri;
/* DOCUMENT yoga_getri
   yoga_getri,d_mat or
   yoga_getri,h_mat, d_mat

   This function computes the inverse of a matrix d_mat (or h_mat) using the LU factorization
   d_mat   : n x n (CArMA obj)
   h_mat   : n x n (yArray)
   
   WARNING : d_mat will be replaced by its inverse
   
   SEE ALSO:
 */

extern yoga_potri;
/* DOCUMENT yoga_potri
   yoga_potri,d_mat

   This function computes the inverse of a real symmetric positive definite matrix mat using the LU factorization
   d_mat   : n x n (CArMA obj)
   
   WARNING : mat will be replaced by its inverse
   
   SEE ALSO:
 */

extern yoga_potri_mgpu;
/* DOCUMENT yoga_potri_mgpu
   yoga_potri_mgpu, ngpus, h_mat, d_mat

   This function computes the inverse of a real symmetric positive definite matrix mat using the LU factorization
   ngpus : number of GPUs
   h_mat   : n x n (yArray)
   d_mat   : n x n (CArMA obj)
   
   WARNING : mat will be replaced by its inverse
   
   SEE ALSO:
 */

extern yoga_svd_host;
/* DOCUMENT yoga_svd_host
   yoga_svd_host,mat,eigenvals,vt,u

   This function computes the svd of matrix mat
   mat   : m x n
   VT    : n x n
   U     : m x m
   eigen : min(m,n)
   
   SEE ALSO:
 */

extern yoga_cula_svd;
/* DOCUMENT yoga_cula_svd
   yoga_cula_svd,mat,eigenvals,vt,u

   This function computes the svd of matrix mat
   mat   : m x n
   VT    : n x n
   U     : m x m
   eigen : min(m,n)
   
   SEE ALSO:
 */

extern yoga_cula_svd_host;
/* DOCUMENT yoga_cula_svd_host
   yoga_cula_svd_host,mat,eigenvals,vt,u

   This function computes the svd of matrix mat
   mat   : m x n
   VT    : n x n
   U     : m x m
   eigen : min(m,n)
   
   SEE ALSO:
 */

/*
__   _____   ____    _      _   _           _      ___  _     _ 
\ \ / / _ \ / ___|  / \    | | | | ___  ___| |_   / _ \| |__ (_)
 \ V / | | | |  _  / _ \   | |_| |/ _ \/ __| __| | | | | '_ \| |
  | || |_| | |_| |/ ___ \  |  _  | (_) \__ \ |_  | |_| | |_) | |
  |_| \___/ \____/_/   \_\ |_| |_|\___/|___/\__|  \___/|_.__// |
                                                           |__/ 
*/

extern yoga_host_obj
/* DOCUMENT yoga_host_obj
   obj = yoga_host_obj("type",[ndims,dims1,dims2,...], pined=1, zero=1)
   or
   obj = yoga_host_obj(0.0f,[ndims,dims1,dims2,...], pined=1, zero=1)
   or
   obj = yoga_host_obj(my_array, pined=1, zero=1)
   or
   obj = yoga_host_obj(other_obj, pined=1, zero=1) 

   This function creates a new Yoga Host object either from a type and
   dimensions list, from a standard Yorick object or from another Yoga Host object
   1rst form: type is : float, double, etc ...
   2nd form: the first argument is a scalarof the same type as the object
   you want to create
   3rd form: my array is a standard Yorick array
   4rth form: other_obj is another yoga object

   Without option, it's doing a stardard memory alloc.
   For speedup transferts, you can pin the memory or doing a "Zero Copy" allocation
   (more details in Jason Sanders, Edward Kandrot - CUDA by Example  An Introduction to General-Purpose GPU Programming - 2010) 
   
   SEE ALSO:
 */

extern yoga_host_getp;
/* DOCUMENT yoga_host_getp
   ptr = yoga_host_getp(yoga_host_obj)

   This function creates a pointer of the yoga_host_obj
   
   SEE ALSO:
 */

extern context_getactivedevice;
extern context_get_maxGflopsDeviceId;

/*
 __  __       _       
|  \/  | __ _(_)_ __  
| |\/| |/ _` | | '_ \ 
| |  | | (_| | | | | |
|_|  |_|\__,_|_|_| |_|
                      
 */

write,"Now loading YoGA version "+yoga_version;
//write, format="yoga_init, %d;\n",useDevice
current_context = yoga_context();
current_context;
yoga_init;
