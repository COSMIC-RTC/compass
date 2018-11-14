#include "naga_obj.hpp"

namespace py = pybind11;
#ifdef CAN_DO_HALF
template void declare_naga_obj<half>(py::module &mod, std::string suffix);

void declare_half_setter_getter(py::module &mod) {
  mod.def("make_naga_obj_half",
          [](carma_context &c,
             const py::array_t<float, py::array::f_style | py::array::forcecast>
                 &data) {
            int ndim = data.ndim() + 1;
            std::vector<long> data_dims(ndim);
            data_dims[0] = data.ndim();
            copy(data.shape(), data.shape() + data.ndim(),
                 begin(data_dims) + 1);
            std::vector<half> tmp(data.size());
            for (int i = 0; i < data.size(); i++) {
              tmp.at(i) = __float2half(data.data()[i]);
            }
            return std::unique_ptr<carma_obj<half>>(new carma_obj<half>(
                &c, data_dims.data(), (const half *)tmp.data()));
          },
          "TODO",  // TODO do the documentation...
          py::arg("context").none(false), py::arg("h_data").none(false));

  mod.def("get_naga_obj_half",
          [](carma_obj<half> &hobj) {
            const long *dims = hobj.getDims();
            std::vector<ssize_t> shape(dims[0]);
            std::vector<ssize_t> strides(dims[0]);
            ssize_t stride = sizeof(float);

            // C-style
            // for (ssize_t dim(dims[0] - 1); dim >= 0; --dim)
            // {
            //   shape[dim] = dims[dim + 1];
            //   strides[dim] = stride;
            //   stride *= shape[dim];
            // }

            // F-style
            for (ssize_t dim(0); dim < dims[0]; ++dim) {
              shape[dim] = dims[dim + 1];
              strides[dim] = stride;
              stride *= shape[dim];
            }

            std::vector<half> tmp(hobj.getNbElem());
            // std::vector<float> tmp2(hobj.getNbElem());
            hobj.device2host(tmp.data());
            py::array_t<float, py::array::f_style | py::array::forcecast> data(
                shape, strides);

            for (int i = 0; i < hobj.getNbElem(); i++) {
              data.mutable_data()[i] = __half2float(tmp.at(i));
            }
            return data;
          },
          "TODO",  // TODO do the documentation...
          py::arg("half_obj").none(false));
}
#endif
