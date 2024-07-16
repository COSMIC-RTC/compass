// This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
//
// COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
// General Public License as published by the Free Software Foundation, either version 3 of the 
// License, or any later version.
//
// COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
// See the GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License along with COMPASS. 
// If not, see <https://www.gnu.org/licenses/>
//
//  Copyright (C) 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>

//! \file      host_obj.cpp
//! \ingroup   libcarma
//! \brief     this file provides pybind wrapper for CarmaHostObj
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#include "host_obj.hpp"
#include <carma.hpp>
#include "type_list.hpp"

using TypeListHostObj =
    GenericTypeList<int32_t, float, double, cuFloatComplex, cuDoubleComplex>;

void declare_carma_host_obj(py::module &mod) {
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
