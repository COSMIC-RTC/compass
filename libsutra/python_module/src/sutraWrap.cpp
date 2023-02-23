// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      sutraWrap.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for sutra
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.3.0
//! \date      2022/01/24

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
void declare_centroider_maskedPix(py::module &);
void declare_controller(py::module &);
void declare_controller_ls(py::module &);
void declare_controller_mv(py::module &);
void declare_controller_geo(py::module &);
void declare_controller_generic(py::module &);
void declare_controller_generic_linear(py::module &);
void declare_controller_cured(py::module &);
void declare_rtc(py::module &);
void declare_gamora(py::module &);
void declare_groot(py::module &);
void declare_coronagraph(py::module &mod);
void declare_perfect_coronagraph(py::module &mod);
void declare_stellar_coronagraph(py::module &mod);

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
  declare_centroider_maskedPix(mod);
  declare_controller(mod);
  declare_controller_ls(mod);
  declare_controller_mv(mod);
  declare_controller_geo(mod);
  declare_controller_generic(mod);
  declare_controller_generic_linear(mod);
  declare_controller_cured(mod);
  declare_rtc(mod);
  declare_gamora(mod);
  declare_groot(mod);
  declare_coronagraph(mod);
  declare_perfect_coronagraph(mod);
  declare_stellar_coronagraph(mod);
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
