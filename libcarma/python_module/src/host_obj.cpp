// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      host_obj.cpp
//! \ingroup   libcarma
//! \brief     this file provides pybind wrapper for CarmaHostObj
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.3.0
//! \date      2022/01/24

#include "host_obj.hpp"
#include <carma.h>
#include "type_list.hpp"

// #ifdef CAN_DO_HALF
// using TypeListHostObj =
//     GenericTypeList<int, float, double, half, cuFloatComplex,
//     cuDoubleComplex>;
// #else
using TypeListHostObj =
    GenericTypeList<int, float, double, cuFloatComplex, cuDoubleComplex>;
// #endif

void declare_carmaWrap_host_obj(py::module &mod) {
  py::enum_<MemAlloc>(mod, "MemAlloc")
      .value("MA_MALLOC", MemAlloc::MA_MALLOC)
      .value("MA_PAGELOCK", MemAlloc::MA_PAGELOCK)
      .value("MA_ZEROCPY", MemAlloc::MA_ZEROCPY)
      .value("MA_PORTABLE", MemAlloc::MA_PORTABLE)
      .value("MA_WRICOMB", MemAlloc::MA_WRICOMB)
      .value("MA_GENEPIN", MemAlloc::MA_GENEPIN)
      .export_values();

  apply<CarmaHostObjInterfacer, TypeListHostObj>(mod);
}
