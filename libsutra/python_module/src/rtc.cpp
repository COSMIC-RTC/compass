#include <wyrm>

#include <sutra_rtc.h>

namespace py = pybind11;

template <typename T>
using F_arrayS = py::array_t<T, py::array::f_style | py::array::forcecast>;
// typedef py::array_t<float, py::array::f_style | py::array::forcecast>
// F_arrayS;

template <typename T>
std::unique_ptr<sutra_rtc<T>> rtc_init() {
  return std::unique_ptr<sutra_rtc<T>>(new sutra_rtc<T>());
};

template <typename T>
void build_cmat_impl(sutra_rtc<T> &sr, int ncontrol, int nfilt, bool filt_tt) {
  sutra_controller_ls<T> *control =
      dynamic_cast<sutra_controller_ls<T> *>(sr.d_control[ncontrol]);
  control->build_cmat(nfilt, filt_tt);
};
#ifdef CAN_DO_HALF
template <>
void build_cmat_impl(sutra_rtc<half> &sr, int ncontrol, int nfilt,
                     bool filt_tt) {
  throw std::runtime_error("Not implemeted");
};
#endif

template <typename T>
void set_gain_impl(sutra_rtc<T> &sr, int ncontrol, T gain) {
  sutra_controller_ls<T> *control =
      dynamic_cast<sutra_controller_ls<T> *>(sr.d_control[ncontrol]);
  control->set_gain(gain);
};
template <>
#ifdef CAN_DO_HALF
void set_gain_impl(sutra_rtc<half> &sr, int ncontrol, half gain) {
  throw std::runtime_error("Not implemeted");
};
#endif

template <typename T>
void set_mgain_impl(sutra_rtc<T> &sr, int ncontrol, F_arrayS<T> data) {
  sutra_controller_ls<T> *control =
      dynamic_cast<sutra_controller_ls<T> *>(sr.d_control[ncontrol]);
  control->set_mgain(data.mutable_data());
};
#ifdef CAN_DO_HALF
template <>
void set_mgain_impl(sutra_rtc<half> &sr, int ncontrol, F_arrayS<half> data) {
  throw std::runtime_error("Not implemeted");
};
#endif

template <typename T>
void svdec_imat_impl(sutra_rtc<T> &sr, int ncontrol) {
  sutra_controller_ls<T> *control =
      dynamic_cast<sutra_controller_ls<T> *>(sr.d_control[ncontrol]);
  control->svdec_imat();
};
#ifdef CAN_DO_HALF
template <>
void svdec_imat_impl(sutra_rtc<half> &sr, int ncontrol) {
  throw std::runtime_error("Not implemeted");
};
#endif

template <typename T>
void rtc_impl(py::module &mod, const char *name) {
  using rtc = sutra_rtc<T>;

  py::class_<rtc>(mod, name)
      .def(py::init(wy::colCast(rtc_init<T>)), " Initialize a void rtc object")

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //

      .def_property_readonly("d_centro",
                             [](rtc &sr) -> vector<sutra_centroider<T> *> & {
                               return sr.d_centro;
                             },
                             "Vector of centroiders")

      .def_property_readonly("d_control",
                             [](rtc &sr) -> vector<sutra_controller<T> *> & {
                               return sr.d_control;
                             },
                             "Vector of controllers")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝
      .def("add_centroider",
           wy::colCast((int (rtc::*)(carma_context *, long, T, T, long, char *,
                                     sutra_wfs *)) &
                       rtc::add_centroider),
           R"pbdoc(
        Add a sutra_centroider object in the RTC

        Parameters
        ------------
        context: (carma_context): carma context
        nvalid:(int): Number of WFS valid ssp
        offset: (float): offset for centroiding computation
        scale: (float): scale factor to get the right unit, ie. arcsec
        device: (int): GPU device index
        typec: (str): Centroider type
        wfs: (sutra_wfs): sutra_wfs handled by the centroider

    )pbdoc",
           py::arg("context"), py::arg("nvalid"), py::arg("offset"),
           py::arg("scale"), py::arg("device"), py::arg("typec"),
           py::arg("wfs"))

      .def("add_centroider",
           wy::colCast(
               (int (rtc::*)(carma_context *, long, T, T, long, char *)) &
               rtc::add_centroider),
           R"pbdoc(
        Add a sutra_centroider object in the RTC

        Parameters
        ------------
        context: (carma_context): carma context
        nvalid:(int): Number of WFS valid ssp
        offset: (float): offset for centroiding computation
        scale: (float): scale factor to get the right unit, ie. arcsec
        device: (int): GPU device index
        typec: (str): Centroider type

    )pbdoc",
           py::arg("context"), py::arg("nvalid"), py::arg("offset"),
           py::arg("scale"), py::arg("device"), py::arg("typec"))

      .def("add_controller", wy::colCast(&rtc::add_controller), R"pbdoc(
        Add a sutra_controller object in the RTC

        Parameters
        ------------
        context: (carma_context): carma context
        nvalid: (int): Number of valid subap.
        nslope: (int): Number of slopes
        nactu:(int): Number of actuators to command
        delay: (float): Loop delay [frames]
        device: (int): GPU device index
        typec: (str): Controller type
        dms: (sutra_dms): sutra_dms object
        idx_dms: (np.array[ndim=1,dtype=np.int64]): index of DM in sutra_dms to command
        ndm: (int): Number of DM to command
        Nphi: (int): Number of pixels in the pupil
        wfs_direction: (bool): Flag for ROKET
        )pbdoc",
           py::arg("context"), py::arg("nvalid"), py::arg("nslope"),
           py::arg("nactu"), py::arg("delay"), py::arg("device"),
           py::arg("typec"), py::arg("dms") = nullptr,
           py::arg("idx_dms") = std::vector<int64_t>(), py::arg("ndm") = 0,
           py::arg("Nphi") = 0, py::arg("wfs_direction") = false)

      .def("do_centroids", (int (rtc::*)(int)) & rtc::do_centroids,
           R"pbdoc(
        Computes the centroids

        Parameters
        ------------
        ncontrol: (int): Index of the controller
    )pbdoc",
           py::arg("ncontrol"))

      .def("do_centroids_geom", &rtc::do_centroids_geom,
           R"pbdoc(
        Computes the centroids geom

        Parameters
        ------------
        ncontrol: (int): Index of the controller
    )pbdoc",
           py::arg("ncontrol"))

      .def("do_centroids_ref", &rtc::do_centroids_ref,
           R"pbdoc(
        Computes the centroids ref

        Parameters
        ------------
        ncontrol: (int): Index of the controller
    )pbdoc",
           py::arg("ncontrol"))

      .def("set_centroids_ref", wy::colCast(&rtc::set_centroids_ref),
           R"pbdoc(
          Set the reference centroids

          Parameters
          ------------
          centroids_ref : (np.array(ndim=1, dtype=np.float32)): ref centroids
     )pbdoc",
           py::arg("centroidsRef"))

      .def("do_control", (int (rtc::*)(int)) & rtc::do_control,
           R"pbdoc(
        Computes the commands

        Parameters
        ------------
        ncontrol: (int): Index of the controller
    )pbdoc",
           py::arg("ncontrol"))

      .def("apply_control", &rtc::apply_control, R"pbdoc(
        Apply the commands on the DM and shape it

        Parameters
        ------------
        ncontrol: (int): Index of the controller
        dms: (sutra_dms): sutra_dms object
        compVoltage: (bool): if True (default), computes delay and perturb voltages. Else, applies just the vector command
    )pbdoc",
           py::arg("ncontrol"), py::arg("dms"), py::arg("compVoltage") = true)

      .def("comp_voltage", &rtc::comp_voltage, R"pbdoc(
        Compute the commands on the DM

        Parameters
        ------------
        ncontrol: (int): Index of the controller
    )pbdoc",
           py::arg("ncontrol"))

      .def("do_imat", wy::colCast(&rtc::do_imat), R"pbdoc(
        Computes interaction matrix

        Parameters
        ------------
        ncontrol: (int): Index of the controller
        dms: (sutra_dms): sutra_dms object
    )pbdoc",
           py::arg("ncontrol"), py::arg("dms"))

      .def("do_imat_basis", wy::colCast(&rtc::do_imat_basis), R"pbdoc(
		Computes a modal interaction matrix

		Parameters
		------------
		ncontrol: (int): Index of the controller
		dms: (sutra_dms): sutra_dms object
        nModes: (int): number of modes in the basis
        m2v: (np.array[ndim=2,dtype=np.float32]): modeToActu matrix
        pushAmpl: (np.array[ndim=1,dtype=np.float32]): pushpull strength in mode units
	)pbdoc",
           py::arg("ncontrol"), py::arg("dms"), py::arg("nModes"),
           py::arg("m2v"), py::arg("pushAmpl"))

      .def("imat_svd", svdec_imat_impl<T>,
           R"pbdoc(
        Computes imat svd

        Parameters
        ------------
        ncontrol : (int): controller index
    )pbdoc",
           py::arg("ncontrol"))

      .def("build_cmat", build_cmat_impl<T>,
           R"pbdoc(
        Computes cmat

        Parameters
        ------------
        ncontrol : (int): controller index
        nfilt: (int): number of modes to filter
        filt_tt: (bool): Flag for TT filter
    )pbdoc",
           py::arg("ncontrol"), py::arg("nfilt"), py::arg("filt_tt") = false)

      .def("do_clipping", &rtc::do_clipping, R"pbdoc(
        Clip the command to apply on the DMs on a sutra_controller object

        Parameters
        ------------
        ncontrol: (int) : controller index
        min: (float) : minimum value for the command
        max: (float) : maximum value for the command
    )pbdoc",
           py::arg("ncontrol"), py::arg("min"), py::arg("max"))

      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
      //
      .def("set_gain", set_gain_impl<T>,
           R"pbdoc(
        Set the loop gain in the controller

        Parameters
        ------------
        ncontrol: (int): controller index
        gain: (float): gain to set
    )pbdoc",
           py::arg("ncontrol"), py::arg("gain"))

      .def("set_mgain", set_mgain_impl<T>,
           R"pbdoc(
        Set the modal gain in the controller

        Parameters
        ------------
        ncontrol: (int): controller index
        ngain: (np.array[ndim=1,dtype=np.float32]): modal gains to set
    )pbdoc",
           py::arg("ncontrol"), py::arg("mgain"));
}

void declare_rtc(py::module &mod) {
#ifdef CAN_DO_HALF
  rtc_impl<half>(mod, "RtcH");
#endif
  rtc_impl<float>(mod, "Rtc");
}
