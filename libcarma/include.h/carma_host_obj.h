// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.
// -----------------------------------------------------------------------------

//! \file      carma_host_obj.h
//! \ingroup   libcarma
//! \class     CarmaHostObj
//! \brief     this class provides wrappers to the generic carma host object
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.4.0
//! \date      2022/01/24

#ifndef _CarmaHostObj_H_
#define _CarmaHostObj_H_

#include <carma.h>
#include <carma_context.h>
#include <carma_streams.h>
#include <carma_utils.h>

#include <iostream>
#include <typeinfo>  // operator typeid

enum MemAlloc {
  MA_MALLOC,
  MA_PAGELOCK,
  MA_ZEROCPY,
  MA_PORTABLE,
  MA_WRICOMB,
  MA_GENEPIN
};

#define MEMORY_ALIGNMENT 4096
#define ALIGN_UP(x, size) (((size_t)x + (size - 1)) & (~(size - 1)))

template <class T_data>
class CarmaObj;

template <class T_data>
class CarmaHostObj {
 protected:
  T_data *h_data;        ///< Input data
  T_data *data_UA;       ///< unpadded input dara for generic pinned mem
  long *dims_data;       ///< dimensions of the array
  int nb_elem;           ///< number of elments in the array
  MemAlloc malloc_type;  ///< type of host alloc
  CarmaStreams *streams;

  void init(const long *dims_data, const T_data *data, MemAlloc malloc_type,
            int nb_streams);

 public:
  CarmaHostObj();
  CarmaHostObj(const long *dims_data);
  CarmaHostObj(const std::vector<long> &dims);
  CarmaHostObj(const long *dims_data, MemAlloc malloc_type);
  CarmaHostObj(const CarmaHostObj<T_data> *obj);
  CarmaHostObj(const CarmaHostObj<T_data> *obj, MemAlloc malloc_type);
  CarmaHostObj(const long *dims_data, const T_data *data);
  CarmaHostObj(const long *dims_data, const T_data *data, MemAlloc malloc_type);
  CarmaHostObj(const long *dims_data, int nb_streams);
  CarmaHostObj(const long *dims_data, MemAlloc malloc_type, int nb_streams);
  CarmaHostObj(const CarmaHostObj<T_data> *obj, int nb_streams);
  CarmaHostObj(const CarmaHostObj<T_data> *obj, MemAlloc malloc_type,
               int nb_streams);
  CarmaHostObj(const long *dims_data, const T_data *data, int nb_streams);
  CarmaHostObj(const long *dims_data, const T_data *data, MemAlloc malloc_type,
               int nb_streams);
  ~CarmaHostObj();

  void get_devpntr(void **pntr_dev);

  int get_nb_streams();
  int add_stream();
  int add_stream(int nb);
  int del_stream();
  int del_stream(int nb);
  cudaStream_t get_cuda_stream(int stream);
  int wait_stream(int stream);
  int wait_all_streams();

  int cpy_obj(CarmaObj<T_data> *carma_obj, cudaMemcpyKind flag);
  int cpy_obj(CarmaObj<T_data> *carma_obj, cudaMemcpyKind flag,
              unsigned int stream);

  /**< General Utilities */
  T_data &operator[](long idx) { return this->h_data[idx]; }
  const T_data &operator[](long idx) const { return this->h_data[idx]; }
  operator T_data *() { return h_data; }
  operator std::string() {
    std::ostringstream stream;
    stream << *this;
    return stream.str();
  }
  // inline char const *c_str() { return std::string(*this).c_str(); }
  T_data *get_data() { return h_data; }
  T_data *get_data_at(int index) { return &h_data[index]; }
  const long *get_dims() { return dims_data; }
  long get_dims(int i) { return dims_data[i]; }
  int get_nb_elements() { return nb_elem; }

  /**< Memory transfer */
  int fill_from(const T_data *data);
  int fill_into(T_data *data);
  int fill(T_data value);

  std::string get_mem_alloc() {
    switch (malloc_type) {
      case MA_MALLOC:
        return "MA_MALLOC";
      case MA_PAGELOCK:
        return "MA_PAGELOCK";
      case MA_ZEROCPY:
        return "MA_ZEROCPY";
      case MA_PORTABLE:
        return "MA_PORTABLE";
      case MA_WRICOMB:
        return "MA_WRICOMB";
      case MA_GENEPIN:
        return "MA_GENEPIN";
      default:
        return "MA_UNKNOWN";
    }
  }
};

template <class T_data>
std::ostream &operator<<(std::ostream &os, CarmaHostObj<T_data> &obj) {
  os << "-----------------------" << std::endl;
  os << "CarmaHostObj<" << typeid(T_data).name() << "> object" << std::endl;
  long ndims = obj.get_dims(0);
  os << "ndims = " << ndims << std::endl;
  for (long dim = 0; dim < ndims; dim++) {
    os << "dim[" << dim << "] = " << obj.get_dims(dim + 1) << std::endl;
  }
  os << "nbElem = " << obj.get_nb_elements() << std::endl;
  os << "sizeof(" << typeid(T_data).name() << ") = " << sizeof(T_data)
     << std::endl;
  os << "-----------------------" << std::endl;
  return os;
}

/*
 extern "C" {

 }
 */

#endif  // _CarmaHostObj_H_
