// This file is part of COMPASS <https://github.com/COSMIC-RTC/compass>
//
// COMPASS is free software: you can redistribute it and/or modify it under the terms of the GNU
// Lesser General Public License as published by the Free Software Foundation, either version 3 of
// the License, or any later version.
//
// COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
// even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License along with COMPASS. If
// not, see <https://www.gnu.org/licenses/>
//
//  Copyright (C) 2011-2024 COSMIC Team <https://github.com/COSMIC-RTC/compass>

//! \file      carma.cpp
//! \ingroup   libcarma
//! \brief     this file provides pybind wrapper for libcarma
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24

#include "declare_name.hpp"

#include <carma.hpp>

namespace py = pybind11;

void declare_carma_host_obj(py::module &);
void declare_carma_obj(py::module &);
void declare_carma_sparse_obj(py::module &);
void declare_carma_context(py::module &);
void declare_carma_timer(py::module &);

// Expose classes and methods to Python
PYBIND11_MODULE(carma, mod) {
  mod.doc() = "";
  auto cDeviceProp = py::class_<cudaDeviceProp>(mod, "cudaDeviceProp");
  declare_carma_context(mod);
  declare_carma_obj(mod);
  declare_carma_host_obj(mod);
  declare_carma_sparse_obj(mod);
  declare_carma_timer(mod);

  // declare_carma_obj<cuDoubleComplex>(mod, "double_complex");

#ifdef VERSION_INFO
#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
  mod.attr("__version__") = TOSTRING(VERSION_INFO);
#else
  mod.attr("__version__") = "dev";
#endif
}
