#include <pybind11/pybind11.h>

namespace py = pybind11;

void declare_shesha_tscreen(py::module &);
void declare_shesha_atmos(py::module &);
void declare_shesha_telescope(py::module &);
void declare_shesha_dms(py::module &);
void declare_shesha_dm(py::module &);
void declare_shesha_source(py::module &);
void declare_shesha_target(py::module &);
void declare_shesha_sensors(py::module &);
void declare_shesha_wfs(py::module &);
void declare_shesha_wfs_sh(py::module &);
// Expose classes and methods to Python
PYBIND11_MODULE(shesha_bind, mod) {
  mod.doc() = "Binding module for libsutra into shesha";

  declare_shesha_tscreen(mod);
  declare_shesha_atmos(mod);
  declare_shesha_telescope(mod);
  declare_shesha_dms(mod);
  declare_shesha_dm(mod);
  declare_shesha_source(mod);
  declare_shesha_target(mod);
  declare_shesha_sensors(mod);
  declare_shesha_wfs(mod);
  declare_shesha_wfs_sh(mod);
}
