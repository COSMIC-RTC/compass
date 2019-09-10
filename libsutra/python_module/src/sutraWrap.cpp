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

//! \file      sutraWrap.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for sutra
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.3.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#include <pybind11/pybind11.h>

namespace py = pybind11;

void declare_tscreen(py::module &);
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
void declare_controller(py::module &);
void declare_controller_ls(py::module &);
void declare_controller_mv(py::module &);
void declare_controller_geo(py::module &);
void declare_controller_generic(py::module &);
void declare_controller_cured(py::module &);
void declare_rtc(py::module &);
void declare_gamora(py::module &);
void declare_groot(py::module &);

#ifdef USE_BRAHMA
void declare_target_brahma(py::module &);
void declare_rtc_brahma(py::module &);
#endif  // USE_BRAHMA

#ifdef USE_CACAO
void declare_rtc_cacao(py::module &);
#endif  // USE_CACAO

// Expose classes and methods to Python
PYBIND11_MODULE(sutraWrap, mod) {
  py::module::import("carmaWrap");
  mod.doc() = "Binding module for libsutra";
  declare_tscreen(mod);
  declare_atmos(mod);
  declare_telescope(mod);
  declare_dms(mod);
  declare_dm(mod);
  declare_source(mod);
  declare_lgs(mod);
  declare_target(mod);
  declare_sensors(mod);
  declare_wfs(mod);
  declare_wfs_sh(mod);
  declare_wfs_pyrhr(mod);
  declare_centroider(mod);
  declare_centroider_tcog(mod);
  declare_centroider_wcog(mod);
  declare_centroider_bpcog(mod);
  declare_centroider_corr(mod);
  declare_centroider_pyr(mod);
  declare_controller(mod);
  declare_controller_ls(mod);
  declare_controller_mv(mod);
  declare_controller_geo(mod);
  declare_controller_generic(mod);
  declare_controller_cured(mod);
  declare_rtc(mod);
  declare_gamora(mod);
  declare_groot(mod);

#ifdef USE_BRAHMA
  declare_target_brahma(mod);
  declare_rtc_brahma(mod);
#endif  // USE_BRAHMA

#ifdef USE_CACAO
  declare_rtc_cacao(mod);
#endif  // USE_CACAO

#ifdef VERSION_INFO
#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
  mod.attr("__version__") = TOSTRING(VERSION_INFO);
#else
  mod.attr("__version__") = "dev";
#endif
}
