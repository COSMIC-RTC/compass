#include <pybind11/pybind11.h>

namespace py = pybind11;

void declare_shesha_tscreen(py::module &);
void declare_shesha_atmos(py::module &);
void declare_shesha_telescope(py::module &);
void declare_shesha_dms(py::module &);
void declare_shesha_dm(py::module &);
void declare_shesha_source(py::module &);
void declare_shesha_lgs(py::module &);
void declare_shesha_target(py::module &);
void declare_shesha_sensors(py::module &);
void declare_shesha_wfs(py::module &);
void declare_shesha_wfs_sh(py::module &);
void declare_shesha_wfs_pyrhr(py::module &);
void declare_shesha_centroider(py::module &);
void declare_shesha_centroider_tcog(py::module &);
void declare_shesha_centroider_wcog(py::module &);
void declare_shesha_centroider_bpcog(py::module &);
void declare_shesha_centroider_corr(py::module &);
void declare_shesha_centroider_pyr(py::module &);
void declare_shesha_controller(py::module &);
void declare_shesha_controller_ls(py::module &);
void declare_shesha_controller_mv(py::module &);
void declare_shesha_controller_geo(py::module &);
void declare_shesha_controller_generic(py::module &);
void declare_shesha_controller_cured(py::module &);
void declare_shesha_rtc(py::module &);
void declare_shesha_gamora(py::module &);

#ifdef USE_BRAHMA
void declare_shesha_target_brahma(py::module &);
void declare_shesha_rtc_brahma(py::module &);
#endif  // USE_BRAHMA

// Expose classes and methods to Python
PYBIND11_MODULE(shesha_bind, mod) {
  mod.doc() = "Binding module for libsutra into shesha";

  declare_shesha_tscreen(mod);
  declare_shesha_atmos(mod);
  declare_shesha_telescope(mod);
  declare_shesha_dms(mod);
  declare_shesha_dm(mod);
  declare_shesha_source(mod);
  declare_shesha_lgs(mod);
  declare_shesha_target(mod);
  declare_shesha_sensors(mod);
  declare_shesha_wfs(mod);
  declare_shesha_wfs_sh(mod);
  declare_shesha_wfs_pyrhr(mod);
  declare_shesha_centroider(mod);
  declare_shesha_centroider_tcog(mod);
  declare_shesha_centroider_wcog(mod);
  declare_shesha_centroider_bpcog(mod);
  declare_shesha_centroider_corr(mod);
  declare_shesha_centroider_pyr(mod);
  declare_shesha_controller(mod);
  declare_shesha_controller_ls(mod);
  declare_shesha_controller_mv(mod);
  declare_shesha_controller_geo(mod);
  declare_shesha_controller_generic(mod);
  declare_shesha_controller_cured(mod);
  declare_shesha_rtc(mod);
  declare_shesha_gamora(mod);

#ifdef USE_BRAHMA
  declare_shesha_target_brahma(mod);
  declare_shesha_rtc_brahma(mod);
#endif  // USE_BRAHMA
}
