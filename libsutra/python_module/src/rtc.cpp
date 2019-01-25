#include <wyrm>

#include <sutra_rtc.h>

namespace py = pybind11;
typedef py::array_t<float, py::array::f_style | py::array::forcecast> F_arrayS;
using rtc = sutra_rtc<float>;
using rtcH = sutra_rtc<half>;

std::unique_ptr<rtc> rtc_init() { return std::unique_ptr<rtc>(new rtc()); };
std::unique_ptr<rtcH> rtcH_init() { return std::unique_ptr<rtcH>(new rtcH()); };

void declare_rtc(py::module &mod) {
  py::class_<rtc>(mod, "Rtc")
      .def(py::init(wy::colCast(rtc_init)), " Initialize a void rtc object")

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //

      .def_property_readonly(
          "d_centro",
          [](rtc &sr) -> vector<sutra_centroider<float> *> & {
            return sr.d_centro;
          },
          "Vector of centroiders")

      .def_property_readonly(
          "d_control",
          [](rtc &sr) -> vector<sutra_controller<float> *> & {
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
           wy::colCast((int (rtc::*)(carma_context *, long, float, float, long,
                                     char *, sutra_wfs *)) &
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
           wy::colCast((int (rtc::*)(carma_context *, long, float, float, long,
                                     char *)) &
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

      .def("imat_svd",
           [](rtc &sr, int ncontrol) {
             sutra_controller_ls<float> *control =
                 dynamic_cast<sutra_controller_ls<float> *>(
                     sr.d_control[ncontrol]);
             control->svdec_imat();
           },
           R"pbdoc(
        Computes imat svd

        Parameters
        ------------
        ncontrol : (int): controller index
    )pbdoc",
           py::arg("ncontrol"))

      .def("build_cmat",
           [](rtc &sr, int ncontrol, int nfilt, bool filt_tt = false) {
             sutra_controller_ls<float> *control =
                 dynamic_cast<sutra_controller_ls<float> *>(
                     sr.d_control[ncontrol]);
             control->build_cmat(nfilt, filt_tt);
           },
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
      .def("set_gain",
           [](rtc &sr, int ncontrol, float gain) {
             sutra_controller_ls<float> *control =
                 dynamic_cast<sutra_controller_ls<float> *>(
                     sr.d_control[ncontrol]);
             control->set_gain(gain);
           },
           R"pbdoc(
        Set the loop gain in the controller

        Parameters
        ------------
        ncontrol: (int): controller index
        gain: (float): gain to set
    )pbdoc",
           py::arg("ncontrol"), py::arg("gain"))

      .def("set_mgain",
           [](rtc &sr, int ncontrol, F_arrayS data) {
             sutra_controller_ls<float> *control =
                 dynamic_cast<sutra_controller_ls<float> *>(
                     sr.d_control[ncontrol]);
             control->set_mgain(data.mutable_data());
           },
           R"pbdoc(
        Set the modal gain in the controller

        Parameters
        ------------
        ncontrol: (int): controller index
        ngain: (np.array[ndim=1,dtype=np.float32]): modal gains to set
    )pbdoc",
           py::arg("ncontrol"), py::arg("mgain"));
};

void declare_rtcH(py::module &mod) {
  py::class_<rtcH>(mod, "RtcH")
      .def(py::init(wy::colCast(rtcH_init)), " Initialize a void rtcH object")

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //

      .def_property_readonly(
          "d_centro",
          [](rtcH &sr) -> vector<sutra_centroider<half> *> & {
            return sr.d_centro;
          },
          "Vector of centroiders")

      .def_property_readonly(
          "d_control",
          [](rtcH &sr) -> vector<sutra_controller<half> *> & {
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
           [](rtcH &sr, carma_context *ctxt, long nvalid, float offset,
              float scale, long device, char *typec) {
             sr.add_centroider(ctxt, nvalid, __float2half(offset),
                               __float2half(scale), device, typec);
           },
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

      .def("add_controller",
           [](rtcH &sr, carma_context *ctxt, int nvalid, int nslope, int nactu,
              float delay, long device, char *typec) {
             sr.add_controller(ctxt, nvalid, nslope, nactu, __float2half(delay),
                               device, typec, nullptr, nullptr, 0, 0, false);
           },
           R"pbdoc(
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
        )pbdoc",
           py::arg("context"), py::arg("nvalid"), py::arg("nslope"),
           py::arg("nactu"), py::arg("delay"), py::arg("device"),
           py::arg("typec"))

      .def("do_centroids", (int (rtcH::*)(int)) & rtcH::do_centroids,
           R"pbdoc(
        Computes the centroids

        Parameters
        ------------
        ncontrol: (int): Index of the controller
    )pbdoc",
           py::arg("ncontrol"))

      .def("do_centroids_geom", &rtcH::do_centroids_geom,
           R"pbdoc(
        Computes the centroids geom

        Parameters
        ------------
        ncontrol: (int): Index of the controller
    )pbdoc",
           py::arg("ncontrol"))

      .def("do_centroids_ref", &rtcH::do_centroids_ref,
           R"pbdoc(
        Computes the centroids ref

        Parameters
        ------------
        ncontrol: (int): Index of the controller
    )pbdoc",
           py::arg("ncontrol"))

      .def("set_centroids_ref", wy::colCast(&rtcH::set_centroids_ref),
           R"pbdoc(
          Set the reference centroids

          Parameters
          ------------
          centroids_ref : (np.array(ndim=1, dtype=np.float16)): ref centroids
     )pbdoc",
           py::arg("centroidsRef"))

      .def("do_control", (int (rtcH::*)(int)) & rtcH::do_control,
           R"pbdoc(
        Computes the commands

        Parameters
        ------------
        ncontrol: (int): Index of the controller
    )pbdoc",
           py::arg("ncontrol"))

      .def("apply_control", &rtcH::apply_control, R"pbdoc(
        Apply the commands on the DM and shape it

        Parameters
        ------------
        ncontrol: (int): Index of the controller
        dms: (sutra_dms): sutra_dms object
        compVoltage: (bool): if True (default), computes delay and perturb voltages. Else, applies just the vector command
    )pbdoc",
           py::arg("ncontrol"), py::arg("dms"), py::arg("compVoltage") = true)

      .def("comp_voltage", &rtcH::comp_voltage, R"pbdoc(
        Compute the commands on the DM

        Parameters
        ------------
        ncontrol: (int): Index of the controller
    )pbdoc",
           py::arg("ncontrol"))

      .def("do_imat", wy::colCast(&rtcH::do_imat), R"pbdoc(
        Computes interaction matrix

        Parameters
        ------------
        ncontrol: (int): Index of the controller
        dms: (sutra_dms): sutra_dms object
    )pbdoc",
           py::arg("ncontrol"), py::arg("dms"))

      .def("do_imat_basis", wy::colCast(&rtcH::do_imat_basis), R"pbdoc(
		Computes a modal interaction matrix

		Parameters
		------------
		ncontrol: (int): Index of the controller
		dms: (sutra_dms): sutra_dms object
        nModes: (int): number of modes in the basis
        m2v: (np.array[ndim=2,dtype=np.float16]): modeToActu matrix
        pushAmpl: (np.array[ndim=1,dtype=np.float16]): pushpull strength in mode units
	)pbdoc",
           py::arg("ncontrol"), py::arg("dms"), py::arg("nModes"),
           py::arg("m2v"), py::arg("pushAmpl"))

      .def("do_clipping",
           [](rtcH &sr, int ncontrol, float min, float max) {
             sr.do_clipping(ncontrol, __float2half(min), __float2half(max));
           },
           R"pbdoc(
        Clip the command to apply on the DMs on a sutra_controller object

        Parameters
        ------------
        ncontrol: (int) : controller index
        min: (float) : minimum value for the command
        max: (float) : maximum value for the command
    )pbdoc",
           py::arg("ncontrol"), py::arg("min"), py::arg("max"));
};
