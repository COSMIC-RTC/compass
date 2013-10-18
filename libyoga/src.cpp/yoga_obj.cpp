#include <yoga_obj.h>

template<class T_data>
yoga_obj<T_data>::yoga_obj(yoga_context *current_context, long *dims_data){
  /** \brief yoga_obj creator.
   * \param current_context : the context in which the yoga_obj is created
   * \param dims_data : the array size Yorick format : [ndims,dims1,dims2,...]
   */
  init(current_context, dims_data, NULL, true, 0);
}

template yObjS::yoga_obj(yoga_context *current_context, long *dims_data);
template yObjD::yoga_obj(yoga_context *current_context, long *dims_data);
template yObjI::yoga_obj(yoga_context *current_context, long *dims_data);
template yObjUI::yoga_obj(yoga_context *current_context, long *dims_data);
template yObjS2::yoga_obj(yoga_context *current_context, long *dims_data);
template yObjD2::yoga_obj(yoga_context *current_context, long *dims_data);

template<class T_data>
yoga_obj<T_data>::yoga_obj(yoga_obj<T_data> *src){
  /** \brief yoga_obj creator.
   * \param src : yoga_obj to copy
   */
  init(src->current_context, src->dims_data, src->d_data, false, src->get_nbStreams());
}

template yObjS::yoga_obj(yObjS *src);
template yObjD::yoga_obj(yObjD *src);
template yObjS2::yoga_obj(yObjS2 *src);
template yObjD2::yoga_obj(yObjD2 *src);

template<class T_data>
yoga_obj<T_data>::yoga_obj(yoga_context *current_context,  yoga_obj<T_data> *src){
  /** \brief yoga_obj creator.
   * \param src : yoga_obj to copy
   */
  init(current_context, src->dims_data, src->d_data, false, 0);
}

template yObjS::yoga_obj(yoga_context *current_context, yObjS *src);
template yObjD::yoga_obj(yoga_context *current_context, yObjD *src);
template yObjS2::yoga_obj(yoga_context *current_context, yObjS2 *src);
template yObjD2::yoga_obj(yoga_context *current_context, yObjD2 *src);

template<class T_data>
yoga_obj<T_data>::yoga_obj(yoga_context *current_context, long *dims_data, T_data *data){
  /** \brief yoga_obj creator.
   * \param dims_data : the array size Yorick format : [ndims,dims1,dims2,...]
   * \param data : the array
   */
  init(current_context, dims_data, data, true, 0);
}

template yObjS::yoga_obj(yoga_context *current_context, long *dims_data, float *data);
template yObjD::yoga_obj(yoga_context *current_context, long *dims_data, double *data);
template yObjI::yoga_obj(yoga_context *current_context, long *dims_data, int *data);
template yObjUI::yoga_obj(yoga_context *current_context, long *dims_data, unsigned int *data);
template yObjS2::yoga_obj(yoga_context *current_context, long *dims_data, float2 *data);
template yObjD2::yoga_obj(yoga_context *current_context, long *dims_data, double2 *data);

template<class T_data>
yoga_obj<T_data>::yoga_obj(yoga_context *current_context, long *dims_data, int nb_streams){
  /** \brief yoga_obj creator.
   * \param dims_data : the array size Yorick format : [ndims,dims1,dims2,...]
   */
  init(current_context, dims_data, NULL, true, nb_streams);
}

template yObjS::yoga_obj(yoga_context *current_context, long *dims_data, int nb_streams);
template yObjD::yoga_obj(yoga_context *current_context, long *dims_data, int nb_streams);
template yObjI::yoga_obj(yoga_context *current_context, long *dims_data, int nb_streams);
template yObjUI::yoga_obj(yoga_context *current_context, long *dims_data, int nb_streams);
template yObjS2::yoga_obj(yoga_context *current_context, long *dims_data, int nb_streams);
template yObjD2::yoga_obj(yoga_context *current_context, long *dims_data, int nb_streams);

template<class T_data>
yoga_obj<T_data>::yoga_obj(yoga_context *current_context, long *dims_data, T_data *data, int nb_streams){
  /** \brief yoga_obj creator.
   * \param dims_data : the array size Yorick format : [ndims,dims1,dims2,...]
   * \param data : the array
   */
  init(current_context, dims_data, data, true, nb_streams);
}

template yObjS::yoga_obj(yoga_context *current_context, long *dims_data, float *data, int nb_streams);
template yObjD::yoga_obj(yoga_context *current_context, long *dims_data, double *data, int nb_streams);
template yObjI::yoga_obj(yoga_context *current_context, long *dims_data, int *data, int nb_streams);
template yObjUI::yoga_obj(yoga_context *current_context, long *dims_data, unsigned int *data, int nb_streams);
template yObjS2::yoga_obj(yoga_context *current_context, long *dims_data, float2 *data, int nb_streams);
template yObjD2::yoga_obj(yoga_context *current_context, long *dims_data, double2 *data, int nb_streams);

template<class T_data>
void yoga_obj<T_data>::init(yoga_context *context, long *dims_data, T_data *data, bool fromHost, int nb_streams){
  this->current_context=context;
  const long size_data = dims_data[0]+1;
  this->dims_data = new long[size_data];
  memcpy(this->dims_data, dims_data, size_data*sizeof(long));

  this->nb_elem = dims_data[1];
  for(int i=2; i<size_data; i++)
    this->nb_elem *= dims_data[i];

  cutilSafeCall(cudaMalloc((void**)&(this->d_data), sizeof(T_data)*this->nb_elem));
  this->device = current_context->get_activeDevice();

#ifdef DEBUG
  printf("YOGA Object created @ %8.8lX on GPU%d\n", (unsigned long)this, current_context->get_activeDevice());
#endif

  if(data==NULL) cutilSafeCall(cudaMemset(this->d_data, 0, sizeof(T_data)*this->nb_elem));

  else if(fromHost) this->host2device(data);
  else this->copyFromDevice(data, this->nb_elem);

  this->plan=0;
  this->gen=0;
  this->d_states = 0;
  this->values = 0;
  this->o_data = 0;
  cutilSafeCall(cudaMalloc((void**)&(this->d_numValid), sizeof(size_t)));

  this->streams = new yoga_streams();
  this->add_stream(nb_streams);

}

/*
template<class T_data>
yoga_obj<T_data>& yoga_obj<T_data>::operator= (const yoga_obj<T_data>& obj){
  if(this->d_data!=0) cutilSafeCall( cudaFree(this->d_data) );
  if(this->dims_data!=0) delete(this->dims_data);
  yoga_obj<T_data> *new_obj = new yoga_obj(yoga_context *current_context, obj);
  return new_obj;
}
*/

template<class T_data> yoga_obj<T_data>::~yoga_obj( ){
  /** \brief yoga_obj destructor.
   */

  //cutilSafeCall( cudaThreadSynchronize() );
  int old_device = current_context->get_activeDevice();
  current_context->set_activeDevice(this->device);

  cutilSafeCall( cudaFree(this->d_data) );
  this->d_data=0;

  delete[](this->dims_data);
  this->dims_data=0;

  cudaFree(this->d_numValid);
  this->d_numValid=0;

  delete this->streams;

#ifdef _USE_CUDPP
  if(this->mScanPlan!=0) {
    cudppDestroyPlan(this->mScanPlan);
  }
#endif

  if (this->values!=0) cutilSafeCall( cudaFree(this->values) );

  //if (this->o_data!=0) cutilSafeCall( cudaFree(this->o_data) );

  /* Destroy the CUFFT plan. */
  if(this->plan != 0) cufftSafeCall( cufftDestroy(this->plan) );

  if (this->gen != 0) this->destroy_prng_host();

  if(this->d_states!=0) cutilSafeCall( cudaFree(this->d_states) );
  this->d_states=0;

#ifdef DEBUG
  printf("YOGA Object deleted @ %8.8lX on GPU%d\n", (unsigned long)this, this->device);
#endif
  current_context->set_activeDevice(old_device);
}

template yObjS::~yoga_obj();
template yObjD::~yoga_obj();
template yObjI::~yoga_obj();
template yObjUI::~yoga_obj();
template yObjD2::~yoga_obj();
template yObjS2::~yoga_obj();

template<class T_data>
int yoga_obj<T_data>::get_nbStreams(){
  /** \brief get the number of streams attached to the host object
   */
	return streams->get_nbStreams();
}

template<class T_data>
int yoga_obj<T_data>::add_stream()
{
	this->streams->add_stream();
	return this->streams->get_nbStreams();
}

template<class T_data>
int yoga_obj<T_data>::add_stream(int nb)
{
	this->streams->add_stream(nb);
	return this->streams->get_nbStreams();
}

template<class T_data>
cudaStream_t yoga_obj<T_data>::get_cudaStream_t(int stream)
{
	return this->streams->get_stream(stream);
}
template cudaStream_t yObjS::get_cudaStream_t(int stream);

template<class T_data>
int yoga_obj<T_data>::del_stream()
{
	this->streams->del_stream();
	return this->streams->get_nbStreams();
}

template<class T_data>
int yoga_obj<T_data>::del_stream(int nb)
{
	this->streams->del_stream(nb);
	return this->streams->get_nbStreams();
}

template<class T_data>
int yoga_obj<T_data>::wait_stream(int stream)
{
	this->streams->wait_stream(stream);
	return EXIT_SUCCESS;
}

template<class T_data>
int yoga_obj<T_data>::wait_all_streams()
{
	this->streams->wait_all_streams();
	return EXIT_SUCCESS;
}
template int yObjS::wait_all_streams();

template<class T_data>
int yoga_obj<T_data>::host2device(T_data *data){
  /** \brief host2device data transfer.
   * \param data : input data
   *
   * this method fills d_input with the imput data
   */
  cutilSafeCall( cudaMemcpy(this->d_data, data, sizeof(T_data)* this->nb_elem, cudaMemcpyHostToDevice) );

  return EXIT_SUCCESS;
}

template int yObjUI::host2device(unsigned int *data);
template int yObjS2::host2device(float2 *data);
template int yObjD2::host2device(double2 *data);

/*
template<class T_data>
T_data* yoga_obj<T_data>::getData(){
  * \brief getData data transfer.
   * \return data : pointer on the data
   *

  return d_data;
}

template float* yObjS::getData();
template double* yObjD::getData();
template unsigned int* yObjUI::getData();
template float2* yObjS2::getData();
template double2* yObjD2::getData();
*/

template<class T_data>
int yoga_obj<T_data>::device2host(T_data *data){
  /** \brief device2host data transfer.
   * \param data : output data
   *
   * this method copies the values in d_output to the output array
   */
  cutilSafeCall( cudaMemcpy(data, this->d_data, sizeof(T_data)* this->nb_elem, cudaMemcpyDeviceToHost) );

  return EXIT_SUCCESS;
}

template int yObjS::device2host(float *data);
template int yObjD::device2host(double *data);
template int yObjI::device2host(int *data);
template int yObjUI::device2host(unsigned int *data);
template int yObjS2::device2host(float2 *data);
template int yObjD2::device2host(double2 *data);

template<class T_data>
int yoga_obj<T_data>::device2hostOpt(T_data *data){
  /** \brief device2host data transfer.
   * \param data : output data
   *
   * this method copies the values in d_output to the output array
   */
	if(this->o_data==0) return EXIT_FAILURE;

  cutilSafeCall( cudaMemcpy(data, this->o_data, sizeof(T_data)* this->nb_elem, cudaMemcpyDeviceToHost) );

  return EXIT_SUCCESS;
}

template int yObjS::device2hostOpt(float *data);
template int yObjD::device2hostOpt(double *data);
template int yObjI::device2hostOpt(int *data);
template int yObjUI::device2hostOpt(unsigned int *data);
template int yObjS2::device2hostOpt(float2 *data);
template int yObjD2::device2hostOpt(double2 *data);

template<class T_data>
int yoga_obj<T_data>::copyToDevice(T_data *data, int nb_elem){
  /** \brief device2host data transfer.
   * \param data : output data
   *
   * this method copies the values in d_output to the output array
   */
  if(nb_elem>this->nb_elem) nb_elem = this->nb_elem;

  cutilSafeCall( cudaMemcpy(data, this->d_data, sizeof(T_data)*nb_elem, cudaMemcpyDeviceToDevice) );

  return EXIT_SUCCESS;
}
template int yObjS::copyToDevice(float *data, int nb_elem);

template<class T_data>
int yoga_obj<T_data>::copyFromDevice(T_data *data, int nb_elem){
  /** \brief device2host data transfer.
   * \param data : output data
   *
   * this method copies the values in d_output to the output array
   */
  if(nb_elem>this->nb_elem) nb_elem = this->nb_elem;

  cutilSafeCall(cudaMemcpy(this->d_data, data, sizeof(T_data)*nb_elem,cudaMemcpyDeviceToDevice) );

  return EXIT_SUCCESS;
}

template<class T_data>
T_data yoga_obj<T_data>::sum(){
  int nBlocks;
  int nThreads;

  sumGetNumBlocksAndThreads(this->nb_elem,0,nBlocks,nThreads);

  reduce<T_data>(this->nb_elem,nThreads,nBlocks,this->d_data,this->d_data);

  cutilCheckMsg("Kernel execution failed");

  // sum partial block sums on GPU
  int s=nBlocks;
  while(s > 1)
    {
      int threads = 0, blocks = 0;
      sumGetNumBlocksAndThreads(s, 0, blocks, threads);

      reduce<T_data>(s, threads, blocks,this->d_data,this->d_data);

      s = (s + (threads*2-1)) / (threads*2);

    }
  T_data* h_odata = new T_data[nBlocks];
  cudaMemcpy(h_odata,this->d_data, sizeof(T_data), cudaMemcpyDeviceToHost);
  return h_odata[0];
}
template float yObjS::sum();
template double yObjD::sum();


template<class T_data>
int yoga_obj<T_data>::transpose(yoga_obj<T_data> *source){

	cutilSafeThreadSync();

  transposeCU(source->d_data,this->d_data,this->dims_data[1] , this->dims_data[2]);

  cutilSafeThreadSync();
  return EXIT_SUCCESS;
}
template int yObjS::transpose(yObjS *data);
template int yObjD::transpose(yObjD *data);


/*
__   __         _      _
\ \ / /__  _ __(_) ___| | __ __      ___ __ __ _ _ __  _ __   ___ _ __ ___
 \ V / _ \| '__| |/ __| |/ / \ \ /\ / / '__/ _` | '_ \| '_ \ / _ \ '__/ __|
  | | (_) | |  | | (__|   <   \ V  V /| | | (_| | |_) | |_) |  __/ |  \__ \
  |_|\___/|_|  |_|\___|_|\_\   \_/\_/ |_|  \__,_| .__/| .__/ \___|_|  |___/
                                                |_|   |_|
 */


int _genericCUS(yObjS *dest,yObjS *source){
  /*! \brief sum wrapper for Yorick (single).
   * \param handle : a pointer to the yoga_obj object we want to sum
   *
   * this function returns the sum of the single-precision array
   */
  return launch_generic1d(source->getData(),dest->getData(),source->getNbElem());
}
int _genericCUD(yObjD *dest,yObjD *source){
  /*! \brief sum wrapper for Yorick (double).
   * \param handle : a pointer to the yoga_obj object we want to sum
   *
   * this function returns the sum of the double-precision array
   */
  return launch_generic1d(source->getData(),dest->getData(),source->getNbElem());
}

int _generic2CUS(yObjS *dest,yObjS *source){
  /*! \brief sum wrapper for Yorick (single).
   * \param handle : a pointer to the yoga_obj object we want to sum
   *
   * this function returns the sum of the single-precision array
   */
  return launch_generic2d(source->getData(),dest->getData(),source->getDims(1),source->getDims(2));
}
int _generic2CUD(yObjD *dest,yObjD *source){
  /*! \brief sum wrapper for Yorick (double).
   * \param handle : a pointer to the yoga_obj object we want to sum
   *
   * this function returns the sum of the double-precision array
   */
  return launch_generic2d(source->getData(),dest->getData(),source->getDims(1),source->getDims(2));
}

