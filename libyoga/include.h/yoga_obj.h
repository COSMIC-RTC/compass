/**
 * \class yoga_obj
 *
 * \ingroup libyoga
 *
 * \brief this class provides wrappers to the generic yoga object
 *
 * \author $Author: dg, as $
 *
 * \version $Revision: 1.0 $
 *
 * \date $Date: 2011/01/28$
 *
 */
#ifndef _YOGA_OBJ_H_
#define _YOGA_OBJ_H_

#ifdef _USE_CUDPP
#include <cudpp.h>
#endif

#include <iostream>
#include <curand_kernel.h>
#include <curand.h>
#include <yoga_utils.h>
#include <yoga_streams.h>
#include <yoga_context.h>

/*
create a memory object
void *memory
int  nb of reference

create a class which contains :
- d_data
- ndims
- dims
- strides
- type

new()

new(existing)


and then
modify yoga_obj so that it is :
an object of the previous class
all the methods of a yoga_obj

 */


#define BLOCK_SZ 16

using namespace std;
enum MemType { MT_DEVICE, MT_DARRAY, MT_HOST, MT_PAGELOCK, MT_ZEROCPY , MT_PORTABLE,
	       MT_WRICOMB, MT_GENEPIN }; //should add texture ?

template<class T_data>
class yoga_data {

protected:
  T_data       *d_data;     ///< Pointer to data
  int          ndims;       ///< Number of dimensions
  int          nb_elem;     ///< Number of elements
  long         *dims_data;  ///< Dimensions
  int          *strides;    ///< Strides for each dimension
  MemType      malloc_type; ///< type of alloc

public:
  T_data  *get_data() { return d_data; }
  int     get_ndims() { return ndims; }
  int     get_nb_elem() { return nb_elem; }
  long    *get_dims_data() { return dims_data; }
  long    get_dims_data(int i) { return dims_data[i]; }
  int     *get_strides() { return strides; }
  int     get_strides(int i) { return strides[i]; }
  MemType get_malloc_type() {return malloc_type;}

};

template<class T_data> class yoga_obj {

 protected:
  T_data *d_data;///< Input data  => change to vector
  T_data *o_data;///< optional data (used for scan / reduction)
  int ndim;
  long *dims_data;///< dimensions of the array
  int nb_elem;///< number of elements in the array
  int device; ///< device where the yoga_obj is allocate
  yoga_context *current_context;

  curandGenerator_t gen;
  curandState *d_states;

  int nThreads;
  int nBlocks;

#ifdef _USE_CUDPP
  CUDPPHandle   mScanPlan;        // CUDPP plan handle for prefix sum
#endif
  bool keysOnly; //< optional flag (used for sort)
  unsigned int *values;///< optional data (used for sort)
  size_t *d_numValid;///< used for compact

  cufftHandle plan;///< FFT plan
  cufftType tPlan;///< FFT plan type

  yoga_streams *streams;

  void init(yoga_context *current_context, long *dims_data, T_data *data, bool fromHost, int nb_streams);

 public:
  yoga_obj(yoga_obj<T_data> *obj);
  yoga_obj(yoga_context *current_context, long *dims_data);
  yoga_obj(yoga_context *current_context, yoga_obj<T_data> *obj);
  yoga_obj(yoga_context *current_context, long *dims_data, T_data *data);
  yoga_obj(yoga_context *current_context, long *dims_data, int nb_streams);
  yoga_obj(yoga_context *current_context, yoga_obj<T_data> *obj, int nb_streams);
  yoga_obj(yoga_context *current_context, long *dims_data, T_data *data, int nb_streams);
  ~yoga_obj();

  int get_nbStreams();
  int add_stream();
  int add_stream(int nb);
  int del_stream();
  int del_stream(int nb);
  cudaStream_t get_cudaStream_t(int stream);
  int wait_stream(int stream);
  int wait_all_streams();

  /**< General Utilities */
  T_data* getData() { return d_data; }
  T_data* getOData() { return o_data; }
  long * getDims() { return dims_data; }
  long getDims(int i) { return dims_data[i]; }
  int getNbElem() { return nb_elem; }
  yoga_context* getContext() { return current_context; }

  bool is_rng_init() { return (gen!=NULL); }


  /**< Memory transfers both ways */
  int host2device(T_data *data);
  int device2host(T_data *data);
  int device2hostOpt(T_data *data);
  int host2deviceVect(T_data *data,int incx, int incy);
  int device2hostVect(T_data *data,int incx,int incy);
  int host2deviceMat(T_data *data,int lda, int ldb);
  int device2hostMat(T_data *data,int lda, int ldb);

  int copyToDevice(T_data *data, int nb_elem);
  int copyFromDevice(T_data *data, int nb_elem);

  cufftHandle* getPlan() {return &plan;};///< FFT plan
  cufftType getTPlan()  {return tPlan;};///< FFT plan type

  unsigned int *getValues() {return values;};///< optional data (used for sort)

  /**< sum */
  T_data sum();

  /**< transpose */
  int transpose(yoga_obj<T_data> *source);
  //yoga_obj<T_data>& operator= (const yoga_obj<T_data>& obj);

  /**< Cublas V2 */
  int imax(int incx);
  int imin(int incx);
  T_data asum(int incx);
  T_data nrm2(int incx);
  T_data dot(yoga_obj<T_data> *source, int incx, int incy);
  void scale(T_data alpha, int incx);
  void swap(yoga_obj<T_data> *source, int incx, int incy);
  void axpy(T_data alpha,yoga_obj<T_data> *source, int incx, int incy);
  void rot(yoga_obj<T_data> *source, int incx, int incy, T_data sc, T_data ss);

  void gemv(char trans, T_data alpha, yoga_obj<T_data> *matA, int lda, yoga_obj<T_data> *vectx, int incx, T_data beta, int incy);

  void ger(T_data alpha, yoga_obj<T_data> *vectx, int incx, yoga_obj<T_data> *vecty, int incy, int lda);

  void gemm(char transa, char transb, T_data alpha, yoga_obj<T_data> *matA, int lda, yoga_obj<T_data> *matB, int ldb, T_data beta, int ldc);

  /**< Curand */
  int init_prng(int device);
  int destroy_prng();
  int prng(T_data *output,char gtype, float alpha, float beta);
  int prng(T_data *output,char gtype, float alpha);
  int prng(char gtype, float alpha,float beta);
  int prng(char gtype, float alpha);
  int prng(char gtype);
  int init_prng_host(int seed);
  int prng_host(char gtype);
  int prng_host(char gtype, T_data alpha);
  int destroy_prng_host();

#ifdef _USE_CUDPP
  /**< cudpp */
  int sort_init(bool keysOnly);
  int sort();
  int config_scan(char *type, int dir, int incl);
  int destroy_scan();
  T_data scan();
  int host2deviceInd(unsigned int *ind);
  int device2hostInd(unsigned int *ind);
  int device2deviceInd(unsigned int *src);
  int compact_init();
  int compact(yoga_obj<T_data> *dest);
  int compact(yoga_obj<T_data> *dest,unsigned int *values);
#else
  /**< fake cudpp */
  int sort_init(bool keysOnly){
		cerr << "!!!!!! CUDPP not compiled !!!!!!"<<endl;
		return EXIT_SUCCESS;
	}
  int sort(){
		cerr << "!!!!!! CUDPP not compiled !!!!!!"<<endl;
		return EXIT_SUCCESS;
	}
  int config_scan(char *type, int dir, int incl){
		cerr << "!!!!!! CUDPP not compiled !!!!!!"<<endl;
		return EXIT_SUCCESS;
	}
  int destroy_scan(){
		cerr << "!!!!!! CUDPP not compiled !!!!!!"<<endl;
		return EXIT_SUCCESS;
	}
  T_data scan(){
		cerr << "!!!!!! CUDPP not compiled !!!!!!"<<endl;
		return EXIT_SUCCESS;
	}
  int host2deviceInd(unsigned int *ind){
		cerr << "!!!!!! CUDPP not compiled !!!!!!"<<endl;
		return EXIT_SUCCESS;
	}
  int device2hostInd(unsigned int *ind){
		cerr << "!!!!!! CUDPP not compiled !!!!!!"<<endl;
		return EXIT_SUCCESS;
	}
  int device2deviceInd(unsigned int *src){
		cerr << "!!!!!! CUDPP not compiled !!!!!!"<<endl;
		return EXIT_SUCCESS;
	}
  int compact_init(){
		cerr << "!!!!!! CUDPP not compiled !!!!!!"<<endl;
		return EXIT_SUCCESS;
	}
  int compact(yoga_obj<T_data> *dest){
		cerr << "!!!!!! CUDPP not compiled !!!!!!"<<endl;
		return EXIT_SUCCESS;
	}
  int compact(yoga_obj<T_data> *dest,unsigned int *values){
		cerr << "!!!!!! CUDPP not compiled !!!!!!"<<endl;
		return EXIT_SUCCESS;
	}
#endif

};

typedef yoga_obj<int> yObjI;
typedef yoga_obj<unsigned int> yObjUI;
typedef yoga_obj<float> yObjS;
typedef yoga_obj<double> yObjD;
typedef yoga_obj<float2> yObjS2;
typedef yoga_obj<double2> yObjD2;
typedef yoga_obj<cuFloatComplex> yObjC;
typedef yoga_obj<cuDoubleComplex> yObjZ;

  // CU functions sum
template <class T_data>
void reduce(int size, int threads, int blocks, T_data *d_idata, T_data *d_odata);

  // CU functions transpose
template<class T_data> int transposeCU(T_data *d_idata,T_data *d_odata,long N1,long N2);

  // CU functions generic
template<class T_data> int launch_generic1d(T_data *d_idata,T_data *d_odata,int N);
template<class T_data> int launch_generic2d(T_data *d_odata,T_data *d_idata,int N1, int N2);

  // CU functions curand
int yoga_prng_init(int *seed,const int nThreads, const int nBlocks, curandState *state);
template <class T> int yoga_prng_cu(T *results,const int nThreads,const int nBlocks,curandState *state,char gtype, int n, float alpha, float beta);

  // CU functions fft
template<class T_in, class T_out> cufftType yoga_select_plan();
template <class T_in, class T_out> void yoga_initfft(long *dims_data,cufftHandle plan,cufftType tPlan);
template <class T_in, class T_out> int yoga_fft(T_in *input, T_out *output, int dir,cufftHandle plan);

  // CU functions generic
template<class T_data> int fillindex(T_data *d_odata,T_data *d_idata,int *indx,int N);
template<class T_data> int fillvalues(T_data *d_odata,unsigned int *indx,int N);
template<class T> int getarray2d(T *d_odata,T *d_idata,int x0, int Ncol,int NC, int N);
template<class T> int fillarray2d(T *d_odata,T *d_idata,int x0, int Ncol,int NC, int N);
template<class T> int fillarray2d2(T *d_odata,T *d_idata,int x0, int Ncol,int NC, int N);
template<class T> int yoga_plus(T *d_odata,T elpha,int N);
template<class T> int yoga_plusai(T *d_odata,T *i_data, int i,int sgn, int N);

  // CU functions fftconv
int fftconv_unpad(float *d_odata,float *d_idata,int fftW, int dataH, int dataW,int N, int n,int nim);
int yoga_initfftconv(yObjS *data_in,yObjS *kernel_in,yObjS *padded_data,yObjC *padded_spectrum,int kernelY, int kernelX);
  // CPP functions fftconv
int yoga_fftconv(yObjS *data_out,yObjS *padded_data,yObjC *padded_spectrum,int kernelY, int kernelX);

#ifdef _USE_MAGMA
  // MAGMA functions
template <class T> int yoga_svd(yoga_obj<T> *imat, yoga_obj<T> *eigenvals, yoga_obj<T> *mod2act, yoga_obj<T> *mes2mod);
#else
template <class T> int yoga_svd(yoga_obj<T> *imat, yoga_obj<T> *eigenvals, yoga_obj<T> *mod2act, yoga_obj<T> *mes2mod){
	cerr << "!!!!!! MAGMA not compiled !!!!!!"<<endl;
	return EXIT_SUCCESS;
}
#endif

// CULA functions
template <class T> int yoga_cula_svd(yoga_obj<T> *imat, yoga_obj<T> *eigenvals, yoga_obj<T> *mod2act, yoga_obj<T> *mes2mod);
//

extern "C" {
  void sumGetNumBlocksAndThreads(int n,int device,int &blocks,int &threads);
  int snapTransformSize(int dataSize);
  bool isPow2(unsigned int x);
  unsigned int nextPow2( unsigned int x );
}

#endif // _YOGA_OBJ_H_
