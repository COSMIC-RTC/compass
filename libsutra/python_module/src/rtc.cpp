// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2019 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.
//  Distributed under GNU - LGPL
//
//  COMPASS is free software: you can redistribute it and/or modify it under the
//  terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or any later version.
//
//  COMPASS: End-to-end AO simulation tool using GPU acceleration
//  The COMPASS platform was designed to meet the need of high-performance for
//  the simulation of AO systems.
//
//  The final product includes a software package for simulating all the
//  critical subcomponents of AO, particularly in the context of the ELT and a
//  real-time core based on several control approaches, with performances
//  consistent with its integration into an instrument. Taking advantage of the
//  specific hardware architecture of the GPU, the COMPASS tool allows to
//  achieve adequate execution speeds to conduct large simulation campaigns
//  called to the ELT.
//
//  The COMPASS platform can be used to carry a wide variety of simulations to
//  both testspecific components of AO of the E-ELT (such as wavefront analysis
//  device with a pyramid or elongated Laser star), and various systems
//  configurations such as multi-conjugate AO.
//
//  COMPASS is distributed in the hope that it will be useful, but WITHOUT ANY
//  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
//  FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
//  details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with COMPASS. If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>.
// -----------------------------------------------------------------------------

//! \file      rtc.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for SutraRtc
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.0.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#include <sutra_rtc.h>

#include <wyrm>

namespace py = pybind11;

template <typename T>
using F_arrayS = py::array_t<T, py::array::f_style | py::array::forcecast>;
// typedef py::array_t<float, py::array::f_style | py::array::forcecast>
// F_arrayS;

template <typename Tin, typename Tcomp, typename Tout>
std::unique_ptr<SutraRtc<Tin, Tcomp, Tout>> rtc_init() {
  return std::unique_ptr<SutraRtc<Tin, Tcomp, Tout>>(
      new SutraRtc<Tin, Tcomp, Tout>());
}

template <typename Tin, typename Tcomp, typename Tout>
void add_controller_impl(SutraRtc<Tin, Tcomp, Tout> &sr, CarmaContext *ctxt,
                         int nvalid, int nslope, int nactu, float delay,
                         long device, std::string typec, int nstates) {
  sr.add_controller(ctxt, nvalid, nslope, nactu, delay, device, typec, nullptr,
                    nullptr, 0, nullptr, 0, 0, false, nstates);
}

template <typename Tin, typename Tcomp, typename Tout>
typename std::enable_if<!std::is_same<Tcomp, half>::value, void>::type
build_cmat_impl(SutraRtc<Tin, Tcomp, Tout> &sr, int ncontrol, int nfilt,
                bool filt_tt) {
  sutra_controller_ls<Tcomp, Tout> *control =
      dynamic_cast<sutra_controller_ls<Tcomp, Tout> *>(sr.d_control[ncontrol]);
  control->build_cmat(nfilt, filt_tt);
}
#ifdef CAN_DO_HALF
template <typename Tin, typename Tcomp, typename Tout>
typename std::enable_if<std::is_same<Tcomp, half>::value, void>::type
build_cmat_impl(SutraRtc<Tin, Tcomp, Tout> &sr, int ncontrol, int nfilt,
                bool filt_tt) {
  throw std::runtime_error("Not implemented");
};
#endif

template <typename Tin, typename Tcomp, typename Tout>
typename std::enable_if<!std::is_same<Tcomp, half>::value, void>::type
set_modal_gains_impl(SutraRtc<Tin, Tcomp, Tout> &sr, int ncontrol,
                     F_arrayS<Tcomp> data) {
  sutra_controller_ls<Tcomp, Tout> *control =
      dynamic_cast<sutra_controller_ls<Tcomp, Tout> *>(sr.d_control[ncontrol]);
  control->set_modal_gains(data.mutable_data());
}
#ifdef CAN_DO_HALF
template <typename Tin, typename Tcomp, typename Tout>
typename std::enable_if<std::is_same<Tcomp, half>::value, void>::type
set_modal_gains_impl(SutraRtc<Tin, Tcomp, Tout> &sr, int ncontrol,
                     F_arrayS<Tcomp> data) {
  throw std::runtime_error("Not implemented");
}
#endif

template <typename Tin, typename Tcomp, typename Tout>
typename std::enable_if<!std::is_same<Tcomp, half>::value, void>::type
svdec_imat_impl(SutraRtc<Tin, Tcomp, Tout> &sr, int ncontrol) {
  sutra_controller_ls<Tcomp, Tout> *control =
      dynamic_cast<sutra_controller_ls<Tcomp, Tout> *>(sr.d_control[ncontrol]);
  control->svdec_imat();
}

#ifdef CAN_DO_HALF
template <typename Tin, typename Tcomp, typename Tout>
typename std::enable_if<std::is_same<Tcomp, half>::value, void>::type
svdec_imat_impl(SutraRtc<Tin, Tcomp, Tout> &sr, int ncontrol) {
  throw std::runtime_error("Not implemented");
}
#endif

template <typename Tin, typename Tcomp, typename Tout>
void rtc_impl(py::module &mod, const char *name) {
  using rtc = SutraRtc<Tin, Tcomp, Tout>;

  py::class_<rtc>(mod, name)
      .def(py::init(wy::colCast(rtc_init<Tin, Tcomp, Tout>)),
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
           wy::colCast((int (rtc::*)(CarmaContext *, long, float, float, bool,
                                     long, std::string, SutraWfs *)) &
                       rtc::add_centroider),

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
           wy::colCast((int (rtc::*)(CarmaContext *, long, float, float, bool,
                                     long, std::string)) &
                       rtc::add_centroider),
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

      .def("add_controller", wy::colCast(&rtc::add_controller), R"pbdoc(
     Add a SutraController object in the RTC

    Args:
        context: (CarmaContext): carma context

        nvalid: (int): Number of valid subap.

        nslope: (int): Number of slopes

        nactu:(int): Number of actuators to command

        delay: (float): Loop delay [frames]

        device: (int): GPU device index

        typec: (str): Controller type

        dms: (SutraDms): SutraDms object

        idx_dms: (np.array[ndim=1,dtype=np.int64]): index of DM in SutraDms to command

        ndm: (int): Number of DM to command

        idx_centro: (np.array[ndim=1,dtype=np.int64]): index of centoiders in sutra_rtc.d_centro to handle

        ncentro: (int): number of centroiders handled

        Nphi: (int): Number of pixels in the pupil

        wfs_direction: (bool): Flag for ROKET

        nstates: (int): (optionnal) number of states
        )pbdoc",
           py::arg("context"), py::arg("nvalid"), py::arg("nslope"),
           py::arg("nactu"), py::arg("delay"), py::arg("device"),
           py::arg("typec"), py::arg("dms") = nullptr,
           py::arg("idx_dms") = std::vector<int64_t>(), py::arg("ndm") = 0,
           py::arg("idx_centro") = std::vector<int64_t>(),
           py::arg("ncentro") = 0, py::arg("Nphi") = 0,
           py::arg("wfs_direction") = false, py::arg("nstates") = 0)

      .def("add_controller", wy::colCast(add_controller_impl<Tin, Tcomp, Tout>),
           R"pbdoc(
     Add a SutraController object in the RTC

    Args:
        context: (CarmaContext): carma context

        nvalid: (int): Number of valid subap.

        nslope: (int): Number of slopes

        nactu:(int): Number of actuators to command

        delay: (float): Loop delay [frames]

        device: (int): GPU device index

        typec: (str): Controller type

        nstates: (int): (optionnal) Number of states
        )pbdoc",
           py::arg("context"), py::arg("nvalid"), py::arg("nslope"),
           py::arg("nactu"), py::arg("delay"), py::arg("device"),
           py::arg("typec"), py::arg("nstates") = 0)

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

      .def("do_calibrate_img", (int (rtc::*)(int)) & rtc::do_calibrate_img,
           R"pbdoc(
     Computes the calibrated image

     Args:
        ncontrol: (int): Index of the controller
    )pbdoc",
           py::arg("ncontrol"))

      .def("do_centroids", (int (rtc::*)(int)) & rtc::do_centroids,
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
    )pbdoc",
           py::arg("ncontrol"))

      .def("do_centroids_ref", &rtc::do_centroids_ref,
           R"pbdoc(
     Computes the centroids ref

     Args:
       ncontrol: (int): Index of the controller
    )pbdoc",
           py::arg("ncontrol"))

      .def("set_centroids_ref", wy::colCast(&rtc::set_centroids_ref),
           R"pbdoc(
     Set the reference centroids

     Args:
          centroids_ref : (np.array(ndim=1, dtype=np.float32)): ref centroids
     )pbdoc",
           py::arg("centroidsRef"))

      .def("do_control", (int (rtc::*)(int)) & rtc::do_control,
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

      .def("do_imat", wy::colCast(&rtc::do_imat), R"pbdoc(
     Computes interaction matrix

     Args:
        ncontrol: (int): Index of the controller

        dms: (SutraDms): SutraDms object
    )pbdoc",
           py::arg("ncontrol"), py::arg("dms"))

      .def("do_imat_basis", wy::colCast(&rtc::do_imat_basis), R"pbdoc(
     Computes a modal interaction matrix

     Args:
          ncontrol: (int): Index of the controller

          dms: (SutraDms): SutraDms object

          nModes: (int): number of modes in the basis

          m2v: (np.array[ndim=2,dtype=np.float32]): modeToActu matrix

          pushAmpl: (np.array[ndim=1,dtype=np.float32]): pushpull strength in mode units
	)pbdoc",
           py::arg("ncontrol"), py::arg("dms"), py::arg("nModes"),
           py::arg("m2v"), py::arg("pushAmpl"))

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

      .def("do_clipping", wy::colCast(&rtc::do_clipping), R"pbdoc(
     Clip the command to apply on the DMs on a SutraController object

     Args:
        ncontrol: (int) : controller index
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
          [](rtc &sr, int ncontrol,
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
