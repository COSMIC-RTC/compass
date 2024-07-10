// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      rtc.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for SutraRtc
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#include "sutraWrapUtils.hpp"

#include <sutra_rtc.hpp>

namespace py = pybind11;

template <typename Tin, typename Tcomp, typename Tout>
std::unique_ptr<SutraRtc<Tin, Tcomp, Tout>> rtc_init() {
  return std::unique_ptr<SutraRtc<Tin, Tcomp, Tout>>(
      new SutraRtc<Tin, Tcomp, Tout>());
}

template <typename Tin, typename Tcomp, typename Tout>
int32_t add_controller(SutraRtc<Tin, Tcomp, Tout> &sr, CarmaContext *context, std::string typec,int64_t device, float delay, int32_t nslope,
               int32_t nactu, int32_t nslope_buffers, int32_t nstates, int32_t nstate_buffers,
               int32_t nmodes, int32_t niir_in, int32_t niir_out, bool polc,
               bool is_modal, SutraDms *dms, ArrayFStyle<int32_t> &idx_dms,
               int32_t ndm,ArrayFStyle<int32_t> &idx_centro, int32_t ncentro, int32_t Nphi,
               bool wfs_direction) {
     return sr.add_controller(context, typec, device, delay, nslope, nactu, nslope_buffers, nstates, nstate_buffers,
                               nmodes,niir_in, niir_out, polc, is_modal,
                               dms, idx_dms.mutable_data(), ndm, idx_centro.mutable_data(), ncentro, Nphi, wfs_direction);
}

template <typename Tin, typename Tcomp, typename Tout>
int32_t set_centroids_ref(SutraRtc<Tin, Tcomp, Tout> &sr, ArrayFStyle<float> &centroids_ref) {
    return sr.set_centroids_ref(centroids_ref.mutable_data());
}

template <typename Tin, typename Tcomp, typename Tout>
int32_t do_imat_basis(SutraRtc<Tin, Tcomp, Tout> &sr, int32_t ncntrl, SutraDms *ydm, int32_t nModes, ArrayFStyle<Tcomp> &m2v,
                  ArrayFStyle<Tcomp> &pushAmpl, int32_t kernconv) {
     return sr.do_imat_basis(
           ncntrl, ydm, nModes, m2v.mutable_data(), pushAmpl.mutable_data(), kernconv);
}

template <typename Tin, typename Tcomp, typename Tout>
void add_controller_impl(SutraRtc<Tin, Tcomp, Tout> &sr, CarmaContext &context,
                    std::string typec, int64_t device, float delay, int32_t nslope, int32_t nactu,
                    int32_t nslope_buffers, int32_t nstates, int32_t nstate_buffers, int32_t nmodes,
                    int32_t niir_in, int32_t niir_out, bool polc,bool is_modal) {
  SutraDms *dms = nullptr;
  int32_t *idx_dms = nullptr;
  int32_t ndm = 0;
  int32_t *idx_centro = nullptr;
  int32_t ncentro = 0;
  int32_t Nphi = 0;
  bool wfs_direction = false;
  sr.add_controller(&context, typec, device, nslope, nslope_buffers, nactu, nstates, nstate_buffers,
                    nmodes,niir_in, niir_out, delay, polc, is_modal,
                    dms, idx_dms, ndm, idx_centro, ncentro, Nphi, wfs_direction);
}

template <typename Tin, typename Tcomp, typename Tout>
typename std::enable_if<!std::is_same<Tcomp, half>::value, void>::type
build_cmat_impl(SutraRtc<Tin, Tcomp, Tout> &sr, int32_t ncontrol, int32_t nfilt,
                bool filt_tt) {
  SutraControllerLs<Tcomp, Tout> *control =
      dynamic_cast<SutraControllerLs<Tcomp, Tout> *>(sr.d_control[ncontrol]);
  control->build_cmat(nfilt, filt_tt);
}
#ifdef CAN_DO_HALF
template <typename Tin, typename Tcomp, typename Tout>
typename std::enable_if<std::is_same<Tcomp, half>::value, void>::type
build_cmat_impl(SutraRtc<Tin, Tcomp, Tout> &sr, int32_t ncontrol, int32_t nfilt,
                bool filt_tt) {
  throw std::runtime_error("Not implemented");
};
#endif

template <typename Tin, typename Tcomp, typename Tout>
typename std::enable_if<!std::is_same<Tcomp, half>::value, void>::type
set_modal_gains_impl(SutraRtc<Tin, Tcomp, Tout> &sr, int32_t ncontrol,
                     ArrayFStyle<Tcomp> data) {
  SutraControllerLs<Tcomp, Tout> *control =
      dynamic_cast<SutraControllerLs<Tcomp, Tout> *>(sr.d_control[ncontrol]);
  control->set_modal_gains(data.mutable_data());
}
#ifdef CAN_DO_HALF
template <typename Tin, typename Tcomp, typename Tout>
typename std::enable_if<std::is_same<Tcomp, half>::value, void>::type
set_modal_gains_impl(SutraRtc<Tin, Tcomp, Tout> &sr, int32_t ncontrol,
                     ArrayFStyle<Tcomp> data) {
  throw std::runtime_error("Not implemented");
}
#endif

template <typename Tin, typename Tcomp, typename Tout>
typename std::enable_if<!std::is_same<Tcomp, half>::value, void>::type
svdec_imat_impl(SutraRtc<Tin, Tcomp, Tout> &sr, int32_t ncontrol) {
  SutraControllerLs<Tcomp, Tout> *control =
      dynamic_cast<SutraControllerLs<Tcomp, Tout> *>(sr.d_control[ncontrol]);
  control->svdec_imat();
}

#ifdef CAN_DO_HALF
template <typename Tin, typename Tcomp, typename Tout>
typename std::enable_if<std::is_same<Tcomp, half>::value, void>::type
svdec_imat_impl(SutraRtc<Tin, Tcomp, Tout> &sr, int32_t ncontrol) {
  throw std::runtime_error("Not implemented");
}
#endif

template <typename Tin, typename Tcomp, typename Tout>
void rtc_impl(py::module &mod, const char *name) {
  using rtc = SutraRtc<Tin, Tcomp, Tout>;

  py::class_<rtc>(mod, name)
      .def(py::init(&rtc_init<Tin, Tcomp, Tout>),
           " Initialize a void rtc object")

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //

      .def_property_readonly(
          "d_centro",
          [](rtc &sr) -> vector<SutraCentroider<Tin, Tcomp> *> & {
            return sr.d_centro;
          },
          "Vector of centroiders")

      .def_property_readonly(
          "d_control",
          [](rtc &sr) -> vector<SutraController<Tcomp, Tout> *> & {
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
           (int32_t (rtc::*)(CarmaContext *, int64_t, float, float, bool,
                                     int64_t, std::string, SutraWfs *)) &
                       rtc::add_centroider,

           R"pbdoc(
     Add a SutraCentroider object in the RTC

    Args:
        context: (CarmaContext): carma context

        nvalid:(int): Number of WFS valid ssp

        offset: (float): offset for centroiding computation

        scale: (float): scale factor to get the right unit, ie. arcsec

        filt_TT: (bool): flag to control TT filtering

        device: (int): GPU device index

        typec: (str): Centroider type

        wfs: (SutraWfs): SutraWfs handled by the centroider

    )pbdoc",
           py::arg("context"), py::arg("nvalid"), py::arg("offset"),
           py::arg("scale"), py::arg("filter_TT"), py::arg("device"),
           py::arg("typec"), py::arg("wfs"))

      .def("add_centroider",
           (int32_t (rtc::*)(CarmaContext *, int64_t, float, float, bool,
                                     int64_t, std::string)) &
                       rtc::add_centroider,
           R"pbdoc(
     Add a SutraCentroider object in the RTC

    Args:
        context: (CarmaContext): carma context

        nvalid:(int): Number of WFS valid ssp

        offset: (float): offset for centroiding computation

        scale: (float): scale factor to get the right unit, ie. arcsec

        filt_TT: (bool): flag to control TT filtering

        device: (int): GPU device index

        typec: (str): Centroider type

    )pbdoc",
           py::arg("context"), py::arg("nvalid"), py::arg("offset"),
           py::arg("scale"), py::arg("filter_TT"), py::arg("device"),
           py::arg("typec"))

      .def("add_controller", &add_controller<Tin, Tcomp, Tout>, R"pbdoc(
     Add a SutraController object in the RTC

    Args:
        context: (CarmaContext): carma context

        typec: (str): Controller type

        device: (int): GPU device index

        delay: (float): Loop delay [frames]

        nslope: (int): Number of slopes

        nactu:(int): Number of actuators to command

    Kwargs:

        nslope_buffers: (int32_t) : Number of historic slopes vectors to use

        nstates       : (int32_t)        : Number of states in state vector

        nstate_buffers: (int32_t)        : Number of historic state vectors to use

        nmodes        : (int32_t)        : Number of modes in mode vector

        niir_in       : (int32_t)        : Number of input mode vectors for iir filter

        niir_out      : (int32_t)        : Number of output mode vectors for iir filter

        polc          : (bool)       : Activate the Pseudo Open Loop Control if available

        is_modal      : (bool)       : Activate projection from modes to actu if available

        dms           : (SutraDms)  : SutraDms object

        idx_dms       : (np.array[ndim=1,dtype=np.int64]) : index of DM in SutraDms to command

        ndm           : (int32_t)        : Number of DM to command

        idx_centro    : (np.array[ndim=1,dtype=np.int64]): (optional) Index of centoiders in sutra_rtc.d_centro to handle

        ncentro       : (int32_t)        : Number of centroiders handled

        Nphi          : (int32_t)        : umber of pixels in the pupil

        wfs_direction : (bool)       : Flag for ROKET

        )pbdoc",
          py::arg("context"), py::arg("typec"), py::arg("device"), py::arg("delay"),
          py::arg("nslope"), py::arg("nactu"),
          py::arg("nslope_buffers")=0, py::arg("nstates")=0, py::arg("nstate_buffers")=0,
          py::arg("nmodes")=0,py::arg("niir_in")=0,py::arg("niir_out")=0,py::arg("polc")=false,
          py::arg("is_modal")=false, py::arg("dms") = nullptr,
          py::arg("idx_dms") = ArrayFStyle<int32_t>(), py::arg("ndm") = 0,
          py::arg("idx_centro") = ArrayFStyle<int32_t>(), py::arg("ncentro") = 0,
          py::arg("Nphi") = 0, py::arg("wfs_direction") = false)

      .def("add_controller", add_controller_impl<Tin, Tcomp, Tout>,
           R"pbdoc(
     Add a SutraController object in the RTC

    Args:
        context: (CarmaContext): carma context

        typec: (str): Controller type

        device: (int): GPU device index

        delay: (float): Loop delay [frames]

        nslope: (int): Number of slopes

        nactu:(int): Number of actuators to command

    Kwargs:
        nslope_buffers: (int32_t) : Number of historic slopes vectors to use

        nstates       : (int32_t)        : Number of states in state vector

        nstate_buffers: (int32_t)        : Number of historic state vectors to use

        nmodes        : (int32_t)        : Number of modes in mode vector

        niir_in       : (int32_t)        : Number of input mode vectors for iir filter

        niir_out      : (int32_t)        : Number of output mode vectors for iir filter

        polc          : (bool)       : Activate the Pseudo Open Loop Control if available

        is_modal      : (bool)       : Activate projection from modes to actu if available

        dms           : (SutraDms)  : SutraDms object

        idx_dms       : (np.array[ndim=1,dtype=np.int64]) : index of DM in SutraDms to command

        ndm           : (int32_t)        : Number of DM to command

        idx_centro    : (np.array[ndim=1,dtype=np.int64]): (optional) Index of centoiders in sutra_rtc.d_centro to handle

        ncentro       : (int32_t)        : Number of centroiders handled

        Nphi          : (int32_t)        : umber of pixels in the pupil

        wfs_direction : (bool)       : Flag for ROKET

        )pbdoc",
          py::arg("context"), py::arg("typec"), py::arg("device"), py::arg("delay"),
          py::arg("nslope"), py::arg("nactu"),
          py::arg("nslope_buffers")=0, py::arg("nstates")=0, py::arg("nstate_buffers")=0,
          py::arg("nmodes")=0,py::arg("niir_in")=0,py::arg("niir_out")=0,py::arg("polc")=false,
          py::arg("is_modal")=false
    )

      .def("remove_controller", &rtc::remove_controller, R"pbdoc(
     Remove the specified controller from the RTC

     Args:
          ncontrol : (int): index of the controller to remove
     )pbdoc",
           py::arg("ncontrol"))

      .def("remove_centroider", &rtc::remove_centroider, R"pbdoc(
     Remove the specified centroider from the RTC

     Args:
          ncentro : (int): index of the centroider to remove
     )pbdoc",
           py::arg("ncentro"))

      .def("do_calibrate_img", (int32_t (rtc::*)(int32_t)) & rtc::do_calibrate_img,
           R"pbdoc(
     Computes the calibrated image

     Args:
        ncontrol: (int): Index of the controller
    )pbdoc",
           py::arg("ncontrol"))

      .def("do_centroids", (int32_t (rtc::*)(int32_t)) & rtc::do_centroids,
           R"pbdoc(
     Computes the centroids

     Args:
        ncontrol: (int): Index of the controller
    )pbdoc",
           py::arg("ncontrol"))

      .def("do_centroids_geom", &rtc::do_centroids_geom,
           R"pbdoc(
     Computes the centroids geom

     Args:
        ncontrol: (int): Index of the controller
        type: (int): Centroiding method of geom wfs
    )pbdoc",
           py::arg("ncontrol"), py::arg("type") = 0)

      .def("do_centroids_ref", &rtc::do_centroids_ref,
           R"pbdoc(
     Computes the centroids ref

     Args:
       ncontrol: (int): Index of the controller
    )pbdoc",
           py::arg("ncontrol"))

      .def("set_centroids_ref", &set_centroids_ref<Tin, Tcomp, Tout>,
           R"pbdoc(
     Set the reference centroids

     Args:
          centroids_ref : (np.array(ndim=1, dtype=np.float32)): ref centroids
     )pbdoc",
           py::arg("centroidsRef"))

      .def("do_control", (int32_t (rtc::*)(int32_t)) & rtc::do_control,
           R"pbdoc(
     Computes the commands

     Args:
       ncontrol: (int): Index of the controller
    )pbdoc",
           py::arg("ncontrol"))

      .def("apply_control", &rtc::apply_control, R"pbdoc(
     Apply the commands on the DM and shape it

     Args:
        ncontrol: (int): Index of the controller

        compVoltage: (bool): if True (default), computes delay and perturb voltages. Else, applies just the vector command
    )pbdoc",
           py::arg("ncontrol"), py::arg("compVoltage") = true)

      .def("comp_voltage", &rtc::comp_voltage, R"pbdoc(
     Compute the commands on the DM

     Args:
        ncontrol: (int): Index of the controller
    )pbdoc",
           py::arg("ncontrol"))

      .def("do_imat", &rtc::do_imat, R"pbdoc(
     Computes interaction matrix

     Args:
        ncontrol: (int): Index of the controller

        dms: (SutraDms): SutraDms object

        kernconv : (bool) : Flag for WFS spot convolution
    )pbdoc",
           py::arg("ncontrol"), py::arg("dms"), py::arg("kernconv"))

      .def("do_imat_basis", &do_imat_basis<Tin, Tcomp, Tout>, R"pbdoc(
     Computes a modal interaction matrix

     Args:
          ncontrol: (int): Index of the controller

          dms: (SutraDms): SutraDms object

          nModes: (int): number of modes in the basis

          m2v: (np.array[ndim=2,dtype=np.float32]): modeToActu matrix

          pushAmpl: (np.array[ndim=1,dtype=np.float32]): pushpull strength in mode units

          kernconv : (bool) : Flag for WFS spot convolution
	)pbdoc",
           py::arg("ncontrol"), py::arg("dms"), py::arg("nModes"),
           py::arg("m2v"), py::arg("pushAmpl"), py::arg("kernconv"))

      .def("imat_svd", svdec_imat_impl<Tin, Tcomp, Tout>,
           R"pbdoc(
     Computes imat svd

     Args:
        ncontrol : (int): controller index
    )pbdoc",
           py::arg("ncontrol"))

      .def("build_cmat", build_cmat_impl<Tin, Tcomp, Tout>,
           R"pbdoc(
     Computes cmat

     Args:
          ncontrol : (int): controller index

          nfilt: (int): number of modes to filter

          filt_tt: (bool): Flag for TT filter
    )pbdoc",
           py::arg("ncontrol"), py::arg("nfilt"), py::arg("filt_tt") = false)

      .def("do_clipping", &rtc::do_clipping, R"pbdoc(
     Clip the command to apply on the DMs on a SutraController object

     Args:
        ncontrol: (int32_t) : controller index
    )pbdoc",
           py::arg("ncontrol"))

      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
      //
      .def(
          "set_gain",
          [](rtc &sr, int32_t ncontrol,
             float gain) { sr.d_control[ncontrol]->set_gain(gain); },
          R"pbdoc(
     Set the loop gain in the controller

     Args:
        ncontrol: (int): controller index

        gain: (float): gain to set
    )pbdoc",
          py::arg("ncontrol"), py::arg("gain"))

      .def("set_modal_gains", set_modal_gains_impl<Tin, Tcomp, Tout>,
           R"pbdoc(
     Set the modal gain in the controller

     Args:
        ncontrol: (int): controller index

        mgain: (np.array[ndim=1,dtype=np.float32]): modal gains to set
    )pbdoc",
           py::arg("ncontrol"), py::arg("mgain"));
}

void declare_rtc(py::module &mod) {
  rtc_impl<float, float, float>(mod, "Rtc_FFF");
  rtc_impl<float, float, uint16_t>(mod, "Rtc_FFU");
  rtc_impl<uint16_t, float, float>(mod, "Rtc_UFF");
  rtc_impl<uint16_t, float, uint16_t>(mod, "Rtc_UFU");
#ifdef CAN_DO_HALF
  rtc_impl<float, half, float>(mod, "Rtc_FHF");
  rtc_impl<float, half, uint16_t>(mod, "Rtc_FHU");
  rtc_impl<uint16_t, half, float>(mod, "Rtc_UHF");
  rtc_impl<uint16_t, half, uint16_t>(mod, "Rtc_UHU");
#endif
}
