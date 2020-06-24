// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2019 COMPASS Team <https://github.com/ANR-COMPASS>
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

//! \file      carma_host_obj.h
//! \ingroup   libcarma
//! \class     carma_host_obj
//! \brief     this class provides wrappers to the generic carma host object
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.4.2
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#ifndef _CARMA_HOST_OBJ_H_
#define _CARMA_HOST_OBJ_H_

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
class carma_obj;

template <class T_data>
class carma_host_obj {
 protected:
  T_data *h_data;       ///< Input data
  T_data *data_UA;      ///< unpadded input dara for generic pinned mem
  long *dims_data;      ///< dimensions of the array
  int nb_elem;          ///< number of elments in the array
  MemAlloc mallocType;  ///< type of host alloc
  carma_streams *streams;

  void init(const long *dims_data, const T_data *data, MemAlloc mallocType,
            int nb_streams);

 public:
  carma_host_obj(const long *dims_data);
  carma_host_obj(const long *dims_data, MemAlloc mallocType);
  carma_host_obj(const carma_host_obj<T_data> *obj);
  carma_host_obj(const carma_host_obj<T_data> *obj, MemAlloc mallocType);
  carma_host_obj(const long *dims_data, const T_data *data);
  carma_host_obj(const long *dims_data, const T_data *data,
                 MemAlloc mallocType);
  carma_host_obj(const long *dims_data, int nb_streams);
  carma_host_obj(const long *dims_data, MemAlloc mallocType, int nb_streams);
  carma_host_obj(const carma_host_obj<T_data> *obj, int nb_streams);
  carma_host_obj(const carma_host_obj<T_data> *obj, MemAlloc mallocType,
                 int nb_streams);
  carma_host_obj(const long *dims_data, const T_data *data, int nb_streams);
  carma_host_obj(const long *dims_data, const T_data *data, MemAlloc mallocType,
                 int nb_streams);
  ~carma_host_obj();

  void get_devpntr(void **pntr_dev);

  int get_nbStreams();
  int add_stream();
  int add_stream(int nb);
  int del_stream();
  int del_stream(int nb);
  cudaStream_t get_cudaStream_t(int stream);
  int wait_stream(int stream);
  int wait_all_streams();

  int cpy_obj(carma_obj<T_data> *caObj, cudaMemcpyKind flag);
  int cpy_obj(carma_obj<T_data> *caObj, cudaMemcpyKind flag,
              unsigned int stream);

  /**< General Utilities */
  operator T_data *() { return h_data; }
  operator std::string() {
    std::ostringstream stream;
    stream << *this;
    return stream.str();
  }
  inline char const *c_str() { return std::string(*this).c_str(); }
  T_data &operator[](int index) { return h_data[index]; }
  T_data *getData() { return h_data; }
  T_data *getDataAt(int index) { return &h_data[index]; }
  const long *getDims() { return dims_data; }
  long getDims(int i) { return dims_data[i]; }
  int getNbElem() { return nb_elem; }

  /**< Memory transfer */
  int fill_from(const T_data *data);
  int fill_into(T_data *data);
  int fill(T_data value);

  std::string getMetAlloc() {
    switch (mallocType) {
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
std::ostream &operator<<(std::ostream &os, carma_host_obj<T_data> &obj) {
  os << "-----------------------" << std::endl;
  os << "carma_host_obj<" << typeid(T_data).name() << "> object" << std::endl;
  long ndims = obj.getDims(0);
  os << "ndims = " << ndims << std::endl;
  for (long dim = 0; dim < ndims; dim++) {
    os << "dim[" << dim << "] = " << obj.getDims(dim + 1) << std::endl;
  }
  os << "nbElem = " << obj.getNbElem() << std::endl;
  os << "sizeof(" << typeid(T_data).name() << ") = " << sizeof(T_data)
     << std::endl;
  os << "-----------------------" << std::endl;
  return os;
}

/*
 extern "C" {

 }
 */

#endif  // _CARMA_HOST_OBJ_H_
