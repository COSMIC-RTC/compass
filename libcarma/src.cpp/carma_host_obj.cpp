#include <carma_host_obj.h>
#include <cuda_runtime.h>
#include <sys/mman.h> // for mmap() / munmap()
template<class T_data>
carma_host_obj<T_data>::carma_host_obj(const long *dims_data) {
  /** \brief carma_host_obj creator.
   * \param dims_data : the array size Yorick format : [ndims,dims1,dims2,...]
   */
  init(dims_data, 0L, MA_MALLOC, 0);
}
template
carma_host_obj<float>::carma_host_obj(const long *dims_data);
template
carma_host_obj<double>::carma_host_obj(const long *dims_data);
template
carma_host_obj<int>::carma_host_obj(const long *dims_data);
template
carma_host_obj<unsigned int>::carma_host_obj(const long *dims_data);
template
carma_host_obj<cuFloatComplex>::carma_host_obj(const long *dims_data);
template
carma_host_obj<cuDoubleComplex>::carma_host_obj(const long *dims_data);

template<class T_data>
carma_host_obj<T_data>::carma_host_obj(const long *dims_data,
    MemAlloc mallocType) {
  /** \brief carma_host_obj creator.
   * \param dims_data : the array size Yorick format : [ndims,dims1,dims2,...]
   * \param mallocType : type of memory allacation : MA_MALLOC, MA_PAGELOCK, MA_ZERO, MA_PORT, MA_WC
   */
  init(dims_data, 0L, mallocType, 0);
}
template
carma_host_obj<float>::carma_host_obj(const long *dims_data,
    MemAlloc mallocType);
template
carma_host_obj<double>::carma_host_obj(const long *dims_data,
    MemAlloc mallocType);
template
carma_host_obj<int>::carma_host_obj(const long *dims_data, MemAlloc mallocType);
template
carma_host_obj<unsigned int>::carma_host_obj(const long *dims_data,
    MemAlloc mallocType);
template
carma_host_obj<cuFloatComplex>::carma_host_obj(const long *dims_data,
    MemAlloc mallocType);
template
carma_host_obj<cuDoubleComplex>::carma_host_obj(const long *dims_data,
    MemAlloc mallocType);

template<class T_data>
carma_host_obj<T_data>::carma_host_obj(carma_host_obj<T_data> *src) {
  /** \brief carma_host_obj creator.
   * \param src : carma_host_obj to copy
   */
  init(src->dims_data, src->h_data, src->mallocType, 0);
}
template
carma_host_obj<float>::carma_host_obj(carma_host_obj<float> *src);
template
carma_host_obj<double>::carma_host_obj(carma_host_obj<double> *src);
template
carma_host_obj<int>::carma_host_obj(carma_host_obj<int> *src);
template
carma_host_obj<unsigned int>::carma_host_obj(carma_host_obj<unsigned int> *src);
template
carma_host_obj<cuFloatComplex>::carma_host_obj(
    carma_host_obj<cuFloatComplex> *src);
template
carma_host_obj<cuDoubleComplex>::carma_host_obj(
    carma_host_obj<cuDoubleComplex> *src);

template<class T_data>
carma_host_obj<T_data>::carma_host_obj(carma_host_obj<T_data> *src,
    MemAlloc mallocType) {
  /** \brief carma_host_obj creator.
   * \param src : carma_host_obj to copy
   * \param mallocType : type of memory allacation : MA_MALLOC, MA_PAGELOCK, MA_ZERO, MA_PORT, MA_WC
   */
  init(src->dims_data, src->h_data, mallocType, 0);
}
template
carma_host_obj<float>::carma_host_obj(carma_host_obj<float> *src,
    MemAlloc mallocType);
template
carma_host_obj<double>::carma_host_obj(carma_host_obj<double> *src,
    MemAlloc mallocType);
template
carma_host_obj<int>::carma_host_obj(carma_host_obj<int> *src,
    MemAlloc mallocType);
template
carma_host_obj<unsigned int>::carma_host_obj(carma_host_obj<unsigned int> *src,
    MemAlloc mallocType);
template
carma_host_obj<cuFloatComplex>::carma_host_obj(
    carma_host_obj<cuFloatComplex> *src, MemAlloc mallocType);
template
carma_host_obj<cuDoubleComplex>::carma_host_obj(
    carma_host_obj<cuDoubleComplex> *src, MemAlloc mallocType);

template<class T_data>
carma_host_obj<T_data>::carma_host_obj(const long *dims_data, T_data *data) {
  /** \brief carma_host_obj creator.
   * \param dims_data : the array size Yorick format : [ndims,dims1,dims2,...]
   * \param data : the array
   */
  init(dims_data, data, MA_MALLOC, 0);
}
template
carma_host_obj<float>::carma_host_obj(const long *dims_data, float *data);
template
carma_host_obj<double>::carma_host_obj(const long *dims_data, double *data);
template
carma_host_obj<int>::carma_host_obj(const long *dims_data, int *data);
template
carma_host_obj<unsigned int>::carma_host_obj(const long *dims_data,
    unsigned int *data);
template
carma_host_obj<cuFloatComplex>::carma_host_obj(const long *dims_data,
    cuFloatComplex *data);
template
carma_host_obj<cuDoubleComplex>::carma_host_obj(const long *dims_data,
    cuDoubleComplex *data);

template<class T_data>
carma_host_obj<T_data>::carma_host_obj(const long *dims_data, T_data *data,
    MemAlloc mallocType) {
  /** \brief carma_host_obj creator.
   * \param dims_data : the array size Yorick format : [ndims,dims1,dims2,...]
   * \param data : the array
   * \param mallocType : type of memory allacation : MA_MALLOC, MA_PAGELOCK, MA_ZERO, MA_PORT, MA_WC
   */
  init(dims_data, data, mallocType, 0);
}
template
carma_host_obj<float>::carma_host_obj(const long *dims_data, float *data,
    MemAlloc mallocType);
template
carma_host_obj<double>::carma_host_obj(const long *dims_data, double *data,
    MemAlloc mallocType);
template
carma_host_obj<int>::carma_host_obj(const long *dims_data, int *data,
    MemAlloc mallocType);
template
carma_host_obj<unsigned int>::carma_host_obj(const long *dims_data,
    unsigned int *data, MemAlloc mallocType);
template
carma_host_obj<cuFloatComplex>::carma_host_obj(const long *dims_data,
    cuFloatComplex *data, MemAlloc mallocType);
template
carma_host_obj<cuDoubleComplex>::carma_host_obj(const long *dims_data,
    cuDoubleComplex *data, MemAlloc mallocType);

template<class T_data>
carma_host_obj<T_data>::carma_host_obj(const long *dims_data, int nb_streams) {
  /** \brief carma_host_obj creator.
   * \param dims_data : the array size Yorick format : [ndims,dims1,dims2,...]
   */
  init(dims_data, 0L, MA_MALLOC, nb_streams);
}
template
carma_host_obj<float>::carma_host_obj(const long *dims_data, int nb_streams);
template
carma_host_obj<double>::carma_host_obj(const long *dims_data, int nb_streams);
template
carma_host_obj<int>::carma_host_obj(const long *dims_data, int nb_streams);
template
carma_host_obj<unsigned int>::carma_host_obj(const long *dims_data,
    int nb_streams);
template
carma_host_obj<cuFloatComplex>::carma_host_obj(const long *dims_data,
    int nb_streams);
template
carma_host_obj<cuDoubleComplex>::carma_host_obj(const long *dims_data,
    int nb_streams);

template<class T_data>
carma_host_obj<T_data>::carma_host_obj(const long *dims_data,
    MemAlloc mallocType, int nb_streams) {
  /** \brief carma_host_obj creator.
   * \param dims_data : the array size Yorick format : [ndims,dims1,dims2,...]
   * \param mallocType : type of memory allacation : MA_MALLOC, MA_PAGELOCK, MA_ZERO, MA_PORT, MA_WC
   */
  init(dims_data, 0L, mallocType, nb_streams);
}
template
carma_host_obj<float>::carma_host_obj(const long *dims_data,
    MemAlloc mallocType, int nb_streams);
template
carma_host_obj<double>::carma_host_obj(const long *dims_data,
    MemAlloc mallocType, int nb_streams);
template
carma_host_obj<int>::carma_host_obj(const long *dims_data, MemAlloc mallocType,
    int nb_streams);
template
carma_host_obj<unsigned int>::carma_host_obj(const long *dims_data,
    MemAlloc mallocType, int nb_streams);
template
carma_host_obj<cuFloatComplex>::carma_host_obj(const long *dims_data,
    MemAlloc mallocType, int nb_streams);
template
carma_host_obj<cuDoubleComplex>::carma_host_obj(const long *dims_data,
    MemAlloc mallocType, int nb_streams);

template<class T_data>
carma_host_obj<T_data>::carma_host_obj(carma_host_obj<T_data> *src,
    int nb_streams) {
  /** \brief carma_host_obj creator.
   * \param src : carma_host_obj to copy
   */
  init(src->dims_data, src->h_data, src->mallocType, nb_streams);
}
template
carma_host_obj<float>::carma_host_obj(carma_host_obj<float> *src,
    int nb_streams);
template
carma_host_obj<double>::carma_host_obj(carma_host_obj<double> *src,
    int nb_streams);
template
carma_host_obj<int>::carma_host_obj(carma_host_obj<int> *src, int nb_streams);
template
carma_host_obj<unsigned int>::carma_host_obj(carma_host_obj<unsigned int> *src,
    int nb_streams);
template
carma_host_obj<cuFloatComplex>::carma_host_obj(
    carma_host_obj<cuFloatComplex> *src, int nb_streams);
template
carma_host_obj<cuDoubleComplex>::carma_host_obj(
    carma_host_obj<cuDoubleComplex> *src, int nb_streams);

template<class T_data>
carma_host_obj<T_data>::carma_host_obj(carma_host_obj<T_data> *src,
    MemAlloc mallocType, int nb_streams) {
  /** \brief carma_host_obj creator.
   * \param src : carma_host_obj to copy
   * \param mallocType : type of memory allacation : MA_MALLOC, MA_PAGELOCK, MA_ZERO, MA_PORT, MA_WC
   */
  init(src->dims_data, src->h_data, mallocType, nb_streams);
}
template
carma_host_obj<float>::carma_host_obj(carma_host_obj<float> *src,
    MemAlloc mallocType, int nb_streams);
template
carma_host_obj<double>::carma_host_obj(carma_host_obj<double> *src,
    MemAlloc mallocType, int nb_streams);
template
carma_host_obj<int>::carma_host_obj(carma_host_obj<int> *src,
    MemAlloc mallocType, int nb_streams);
template
carma_host_obj<unsigned int>::carma_host_obj(carma_host_obj<unsigned int> *src,
    MemAlloc mallocType, int nb_streams);
template
carma_host_obj<cuFloatComplex>::carma_host_obj(
    carma_host_obj<cuFloatComplex> *src, MemAlloc mallocType, int nb_streams);
template
carma_host_obj<cuDoubleComplex>::carma_host_obj(
    carma_host_obj<cuDoubleComplex> *src, MemAlloc mallocType, int nb_streams);

template<class T_data>
carma_host_obj<T_data>::carma_host_obj(const long *dims_data, T_data *data,
    int nb_streams) {
  /** \brief carma_host_obj creator.
   * \param dims_data : the array size Yorick format : [ndims,dims1,dims2,...]
   * \param data : the array
   */
  init(dims_data, data, MA_MALLOC, nb_streams);
}
template
carma_host_obj<float>::carma_host_obj(const long *dims_data, float *data,
    int nb_streams);
template
carma_host_obj<double>::carma_host_obj(const long *dims_data, double *data,
    int nb_streams);
template
carma_host_obj<int>::carma_host_obj(const long *dims_data, int *data,
    int nb_streams);
template
carma_host_obj<unsigned int>::carma_host_obj(const long *dims_data,
    unsigned int *data, int nb_streams);
template
carma_host_obj<cuFloatComplex>::carma_host_obj(const long *dims_data,
    cuFloatComplex *data, int nb_streams);
template
carma_host_obj<cuDoubleComplex>::carma_host_obj(const long *dims_data,
    cuDoubleComplex *data, int nb_streams);

template<class T_data>
carma_host_obj<T_data>::carma_host_obj(const long *dims_data, T_data *data,
    MemAlloc mallocType, int nb_streams) {
  /** \brief carma_host_obj creator.
   * \param dims_data : the array size Yorick format : [ndims,dims1,dims2,...]
   * \param data : the array
   * \param mallocType : type of memory allacation : MA_MALLOC, MA_PAGELOCK, MA_ZERO, MA_PORT, MA_WC
   */
  init(dims_data, data, mallocType, nb_streams);
}
template
carma_host_obj<float>::carma_host_obj(const long *dims_data, float *data,
    MemAlloc mallocType, int nb_streams);
template
carma_host_obj<double>::carma_host_obj(const long *dims_data, double *data,
    MemAlloc mallocType, int nb_streams);
template
carma_host_obj<int>::carma_host_obj(const long *dims_data, int *data,
    MemAlloc mallocType, int nb_streams);
template
carma_host_obj<unsigned int>::carma_host_obj(const long *dims_data,
    unsigned int *data, MemAlloc mallocType, int nb_streams);
template
carma_host_obj<cuFloatComplex>::carma_host_obj(const long *dims_data,
    cuFloatComplex *data, MemAlloc mallocType, int nb_streams);
template
carma_host_obj<cuDoubleComplex>::carma_host_obj(const long *dims_data,
    cuDoubleComplex *data, MemAlloc mallocType, int nb_streams);

template<class T_data>
void carma_host_obj<T_data>::init(const long *dims_data, T_data *data,
    MemAlloc mallocType, int nb_streams) {
  const long size_data = dims_data[0] + 1;
  this->dims_data = new long[size_data];
  copy(dims_data, dims_data + size_data, this->dims_data);

  this->nb_elem = dims_data[1];
  for (int i = 2; i < size_data; i++)
    this->nb_elem *= dims_data[i];

  this->mallocType = mallocType;
  if (mallocType == MA_MALLOC) {
    h_data = new T_data[this->nb_elem];
  } else if (mallocType == MA_PAGELOCK) {
    cutilSafeCall(
        cudaHostAlloc((void** )&(this->h_data), sizeof(T_data) * this->nb_elem, cudaHostAllocDefault));
  } else if (mallocType == MA_ZEROCPY) {
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    if (!prop.canMapHostMemory) {
      cerr << "Can't map host memory\n";
      throw "Can't map host memory\n";
    }
    cudaSetDeviceFlags(cudaDeviceMapHost);
    cutilSafeCall(
        cudaHostAlloc((void** )&(this->h_data), sizeof(T_data) * this->nb_elem, cudaHostAllocWriteCombined | cudaHostAllocMapped));
  } else if (mallocType == MA_PORTABLE) {
    cutilSafeCall(
        cudaHostAlloc((void** )&(this->h_data), sizeof(T_data) * this->nb_elem, cudaHostAllocWriteCombined | cudaHostAllocMapped | cudaHostAllocPortable));
  } else if (mallocType == MA_WRICOMB) {
    cutilSafeCall(
        cudaHostAlloc((void** )&(this->h_data), sizeof(T_data) * this->nb_elem, cudaHostAllocWriteCombined));
  } else if (mallocType == MA_GENEPIN) {
    cudaSetDeviceFlags(cudaDeviceBlockingSync | cudaDeviceMapHost);
    this->data_UA = (T_data *) mmap(NULL,
        (sizeof(T_data) * this->nb_elem + MEMORY_ALIGNMENT),
        PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANON, -1, 0);
    this->h_data = (T_data *) ALIGN_UP( data_UA, MEMORY_ALIGNMENT );
    cutilSafeCall(
        cudaHostRegister(h_data, sizeof(T_data) * this->nb_elem, cudaHostRegisterMapped));
  } else
    throw "Error : type of malloc unknown";

  streams = new carma_streams();
  if (mallocType != MA_MALLOC) {
    this->add_stream(nb_streams);
  }

#if DEBUG  
  printf("CARMA Host Object created @ %8.8lX with %s memory\n", (unsigned long)this, getMetAlloc().c_str());
#endif
  if (data == 0L) {
    if (mallocType == MA_MALLOC) {
      memset(this->h_data, 0, sizeof(T_data) * this->nb_elem);
    } else {
      //cutilSafeCall(cudaMemset((T_data *)this->h_data, 0, sizeof(T_data)*this->nb_elem));
      memset(this->h_data, 0, sizeof(T_data) * this->nb_elem);
    }
  } else
    this->fill_from(data);
}

template<class T_data>
void carma_host_obj<T_data>::get_devpntr(void **pntr_dev) {
  /** \brief retreive device pointer from a host object
   */
  //cutilSafeCall(cudaHostGetDevicePointer((void **)&d_a, (void *)(this->d_data), this->mallocType));
  cutilSafeCall(cudaHostGetDevicePointer(pntr_dev, (void * )(this->h_data), 0));
}

/*
 template<class T_data>
 carma_host_obj<T_data>& carma_host_obj<T_data>::operator= (const carma_host_obj<T_data>& obj){
 if(this->d_data!=0L) cutilSafeCall( cudaFree(this->d_data) );
 if(this->dims_data!=0L) delete(this->dims_data);
 carma_host_obj<T_data> *new_obj = new carma_host_obj(obj);
 return new_obj;
 }
 */

template<class T_data>
carma_host_obj<T_data>::~carma_host_obj() {
  /** \brief carma_host_obj destructor.
   */
  //cutilSafeCall( cudaThreadSynchronize() );
  if (this->h_data != 0L) {
    if (mallocType == MA_MALLOC)
      delete[] this->h_data;
    else if (mallocType == MA_GENEPIN) {
      cutilSafeCall(cudaHostUnregister(this->h_data));
      munmap(this->data_UA, sizeof(T_data) * this->nb_elem);
    } else
      cutilSafeCall(cudaFreeHost(this->h_data));
    this->h_data = 0L;
  }

  delete this->streams;

#if DEBUG
  printf("CARMA Host Object deleted @ %8.8lX\n", (unsigned long)this);
#endif
}
template
carma_host_obj<float>::~carma_host_obj();
template
carma_host_obj<double>::~carma_host_obj();
template
carma_host_obj<int>::~carma_host_obj();
template
carma_host_obj<unsigned int>::~carma_host_obj();
template
carma_host_obj<cuFloatComplex>::~carma_host_obj();
template
carma_host_obj<cuDoubleComplex>::~carma_host_obj();

template<class T_data>
int carma_host_obj<T_data>::get_nbStreams() {
  /** \brief get the number of streams attached to the host object
   */
  return streams->get_nbStreams();
}

template<class T_data>
int carma_host_obj<T_data>::add_stream() {
  if (this->mallocType != MA_MALLOC)
    this->streams->add_stream();
  else
    cerr << "malloc Carma Host Object can not be streamed\n";
  return this->streams->get_nbStreams();
}

template<class T_data>
int carma_host_obj<T_data>::add_stream(int nb) {
  if (this->mallocType != MA_MALLOC)
    this->streams->add_stream(nb);
  else
    cerr << "malloc Carma Host Object can not be streamed\n";
  return this->streams->get_nbStreams();
}

template<class T_data>
cudaStream_t carma_host_obj<T_data>::get_cudaStream_t(int stream) {
  return this->streams->get_stream(stream);
}

template<class T_data>
int carma_host_obj<T_data>::del_stream() {
  this->streams->del_stream();
  return this->streams->get_nbStreams();
}

template<class T_data>
int carma_host_obj<T_data>::del_stream(int nb) {
  this->streams->del_stream(nb);
  return this->streams->get_nbStreams();
}

template<class T_data>
int carma_host_obj<T_data>::wait_stream(int stream) {
  this->streams->wait_stream(stream);
  return EXIT_SUCCESS;
}

template<class T_data>
int carma_host_obj<T_data>::wait_all_streams() {
  this->streams->wait_all_streams();
  return EXIT_SUCCESS;
}

/*
 template<class T_data>
 T_data* carma_host_obj<T_data>::getData(){
 * \brief getData data transfer.
 * \return data : pointer on the data
 *

 return h_data;
 }
 */

template<class T_data>
int carma_host_obj<T_data>::fill_from(T_data *data) {
  if (get_nbStreams() > 1) {
    int nstreams = get_nbStreams();
    for (int i = 0; i < nstreams; i++) {
      cutilSafeCall(
          cudaMemcpyAsync(&(h_data[i * nb_elem / nstreams]),
              &(data[i * nb_elem / nstreams]),
              sizeof(T_data) * nb_elem / nstreams, cudaMemcpyHostToHost,
              get_cudaStream_t(i)));
    }
  } else if (get_nbStreams() == 1) {
    cutilSafeCall(
        cudaMemcpyAsync(h_data, data, nb_elem * sizeof(T_data),
            cudaMemcpyHostToHost, get_cudaStream_t(0)));
  } else {
    if (mallocType == MA_MALLOC) {
      copy(data, data + nb_elem, h_data);
    } else {
      cutilSafeCall(
          cudaMemcpy(h_data, data, nb_elem * sizeof(T_data),
              cudaMemcpyHostToHost));
    }
  }
  return EXIT_SUCCESS;
}
template int
carma_host_obj<float>::fill_from(float *data);
template int
carma_host_obj<double>::fill_from(double *data);
template int
carma_host_obj<int>::fill_from(int *data);
template int
carma_host_obj<unsigned int>::fill_from(unsigned int *data);
template int
carma_host_obj<cuFloatComplex>::fill_from(cuFloatComplex *data);
template int
carma_host_obj<cuDoubleComplex>::fill_from(cuDoubleComplex *data);

template<class T_data>
int carma_host_obj<T_data>::fill_into(T_data *data) {
  if (get_nbStreams() > 1) {
    int nstreams = get_nbStreams();
    for (int i = 0; i < nstreams; i++) {
      cutilSafeCall(
          cudaMemcpyAsync(&(data[i * nb_elem / nstreams]),
              &(h_data[i * nb_elem / nstreams]),
              sizeof(T_data) * nb_elem / nstreams, cudaMemcpyHostToHost,
              get_cudaStream_t(i)));
    }
  } else if (get_nbStreams() == 1) {
    cutilSafeCall(
        cudaMemcpyAsync(data, h_data, nb_elem * sizeof(T_data),
            cudaMemcpyHostToHost, get_cudaStream_t(0)));
  } else {
    if (mallocType == MA_MALLOC) {
      copy(h_data, h_data + nb_elem, data);
    } else {
      cutilSafeCall(
          cudaMemcpy(data, h_data, nb_elem * sizeof(T_data),
              cudaMemcpyHostToHost));
    }
  }
  return EXIT_SUCCESS;
}
template int
carma_host_obj<float>::fill_into(float *data);
template int
carma_host_obj<double>::fill_into(double *data);
template int
carma_host_obj<int>::fill_into(int *data);
template int
carma_host_obj<unsigned int>::fill_into(unsigned int *data);
template int
carma_host_obj<cuFloatComplex>::fill_into(cuFloatComplex *data);
template int
carma_host_obj<cuDoubleComplex>::fill_into(cuDoubleComplex *data);

template<class T_data>
int carma_host_obj<T_data>::cpy_obj(carma_obj<T_data> *caObj,
    cudaMemcpyKind flag) {
  if (nb_elem != caObj->getNbElem()) {
    cerr << "***** ERROR (" << __FILE__ << "@" << __LINE__
        << ") : two objects do not have the same size *****\n";
    return EXIT_FAILURE;
  }

  T_data *data_src, *data_dst;
  if (flag == cudaMemcpyHostToDevice) {
    data_dst = *caObj;
    data_src = h_data;
  } else if (flag == cudaMemcpyDeviceToHost) {
    data_dst = h_data;
    data_src = *caObj;
  } else {
    cerr << "***** ERROR (" << __FILE__ << "@" << __LINE__
        << ") : wrong flag *****\n";
    return EXIT_FAILURE;
  }

  if (get_nbStreams() > 1) {
    int nstreams = get_nbStreams();
    for (int i = 0; i < nstreams; i++) {
      cutilSafeCall(
          cudaMemcpyAsync(&(data_dst[i * nb_elem / nstreams]),
              &(data_src[i * nb_elem / nstreams]),
              sizeof(T_data) * nb_elem / nstreams, flag, get_cudaStream_t(i)));
    }
  } else if (get_nbStreams() == 1) {
    cutilSafeCall(
        cudaMemcpyAsync(data_dst, data_src, sizeof(T_data) * nb_elem, flag,
            get_cudaStream_t(0)));
  } else {
    cutilSafeCall(
        cudaMemcpy(data_dst, data_src, nb_elem * sizeof(T_data), flag));
  }
  return EXIT_SUCCESS;
}
template int
carma_host_obj<float>::cpy_obj(carma_obj<float> *data, cudaMemcpyKind flag);
template int
carma_host_obj<double>::cpy_obj(carma_obj<double> *data, cudaMemcpyKind flag);
template int
carma_host_obj<int>::cpy_obj(carma_obj<int> *data, cudaMemcpyKind flag);
template int
carma_host_obj<unsigned int>::cpy_obj(carma_obj<unsigned int> *data,
    cudaMemcpyKind flag);
template int
carma_host_obj<cuFloatComplex>::cpy_obj(carma_obj<cuFloatComplex> *data,
    cudaMemcpyKind flag);
template int
carma_host_obj<cuDoubleComplex>::cpy_obj(carma_obj<cuDoubleComplex> *data,
    cudaMemcpyKind flag);

template<class T_data>
int carma_host_obj<T_data>::cpy_obj(carma_obj<T_data> *caObj,
    cudaMemcpyKind flag, unsigned int stream) {
  unsigned int nbStreams = this->streams->get_nbStreams();
  if (stream >= nbStreams)
    return EXIT_FAILURE;

  carma_streams streams_tmp = *(this->streams);
  if (flag == cudaMemcpyHostToDevice) {
    cutilSafeCall(
        cudaMemcpyAsync(*caObj, h_data, caObj->getNbElem() * sizeof(T_data),
            flag, get_cudaStream_t(stream)));
  } else if (flag == cudaMemcpyDeviceToHost) {
    cutilSafeCall(
        cudaMemcpyAsync(h_data, *caObj, caObj->getNbElem() * sizeof(T_data),
            flag, get_cudaStream_t(stream)));
  } else
    return EXIT_FAILURE;
  return EXIT_SUCCESS;
}
template int
carma_host_obj<float>::cpy_obj(carma_obj<float> *data, cudaMemcpyKind flag,
    unsigned int stream);
template int
carma_host_obj<double>::cpy_obj(carma_obj<double> *data, cudaMemcpyKind flag,
    unsigned int stream);
template int
carma_host_obj<int>::cpy_obj(carma_obj<int> *data, cudaMemcpyKind flag,
    unsigned int stream);
template int
carma_host_obj<unsigned int>::cpy_obj(carma_obj<unsigned int> *data,
    cudaMemcpyKind flag, unsigned int stream);
template int
carma_host_obj<cuFloatComplex>::cpy_obj(carma_obj<cuFloatComplex> *data,
    cudaMemcpyKind flag, unsigned int stream);
template int
carma_host_obj<cuDoubleComplex>::cpy_obj(carma_obj<cuDoubleComplex> *data,
    cudaMemcpyKind flag, unsigned int stream);

