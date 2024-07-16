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

//! \file      sutra.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for sutra
//! \author    COSMIC Team <https://github.com/COSMIC-RTC/compass>
//! \date      2022/01/24
#define PYBIND11_DETAILED_ERROR_MESSAGES

#include "sutraUtils.hpp"

#include <carma.hpp>

namespace py = pybind11;

void declare_turbu_screen(py::module &);
void declare_atmos(py::module &);
void declare_telescope(py::module &);
void declare_dms(py::module &);
void declare_dm(py::module &);
void declare_source(py::module &);
void declare_lgs(py::module &);
void declare_target(py::module &);
void declare_sensors(py::module &);
void declare_wfs(py::module &);
void declare_wfs_sh(py::module &);
void declare_wfs_pyrhr(py::module &);
void declare_centroider(py::module &);
void declare_centroider_tcog(py::module &);
void declare_centroider_wcog(py::module &);
void declare_centroider_bpcog(py::module &);
void declare_centroider_corr(py::module &);
void declare_centroider_pyr(py::module &);
void declare_centroider_maskedPix(py::module &);
void declare_controller(py::module &);
void declare_controller_ls(py::module &);
void declare_controller_mv(py::module &);
void declare_controller_geo(py::module &);
void declare_controller_generic(py::module &);
void declare_controller_generic_linear(py::module &);
void declare_rtc(py::module &);
void declare_gamora(py::module &);
void declare_groot(py::module &);
void declare_coronagraph(py::module &mod);
void declare_perfect_coronagraph(py::module &mod);
void declare_stellar_coronagraph(py::module &mod);

// Expose classes and methods to Python
PYBIND11_MODULE(sutra, mod) {
  mod.doc() = "Binding module for libsutra";
  declare_turbu_screen(mod);
  declare_atmos(mod);
  declare_telescope(mod);
  declare_dm(mod);
  declare_dms(mod);
  declare_lgs(mod);
  declare_source(mod);
  declare_target(mod);
  declare_wfs(mod);
  declare_wfs_sh(mod);
  declare_wfs_pyrhr(mod);
  declare_sensors(mod);
  declare_centroider(mod);
  declare_centroider_tcog(mod);
  declare_centroider_wcog(mod);
  declare_centroider_bpcog(mod);
  declare_centroider_corr(mod);
  declare_centroider_pyr(mod);
  declare_centroider_maskedPix(mod);
  declare_controller(mod);
  declare_controller_ls(mod);
  declare_controller_mv(mod);
  declare_controller_geo(mod);
  declare_controller_generic(mod);
  declare_controller_generic_linear(mod);
  declare_rtc(mod);
  declare_gamora(mod);
  declare_groot(mod);
  declare_coronagraph(mod);
  declare_perfect_coronagraph(mod);
  declare_stellar_coronagraph(mod);

#ifdef VERSION_INFO
#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
  mod.attr("__version__") = TOSTRING(VERSION_INFO);
#else
  mod.attr("__version__") = "dev";
#endif
}
