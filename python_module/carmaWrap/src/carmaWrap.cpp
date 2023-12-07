// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      carmaWrap.cpp
//! \ingroup   libcarma
//! \brief     this file provides pybind wrapper for libcarma
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#include "declare_name.hpp"

#include <carma.h>

namespace py = pybind11;

void declare_carmaWrap_host_obj(py::module &);
void declare_carmaWrap_obj(py::module &);
void declare_carmaWrap_sparse_obj(py::module &);
void declare_carmaWrap_context(py::module &);
void declare_carmaWrap_timer(py::module &);

// Expose classes and methods to Python
PYBIND11_MODULE(carmaWrap, mod) {
  mod.doc() = "";
  auto cDeviceProp = py::class_<cudaDeviceProp>(mod, "cudaDeviceProp");
  declare_carmaWrap_context(mod);
  declare_carmaWrap_obj(mod);
  declare_carmaWrap_host_obj(mod);
  declare_carmaWrap_sparse_obj(mod);
  declare_carmaWrap_timer(mod);

  // declare_carmaWrap_obj<cuDoubleComplex>(mod, "double_complex");

#ifdef VERSION_INFO
#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
  mod.attr("__version__") = TOSTRING(VERSION_INFO);
#else
  mod.attr("__version__") = "dev";
#endif
}
