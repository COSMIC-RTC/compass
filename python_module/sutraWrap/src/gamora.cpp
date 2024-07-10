// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2023 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      gamora.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for SutraGamora
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.5.0
//! \date      2022/01/24

#include "sutraWrapUtils.hpp"

#include <sutra_gamora.hpp>

namespace py = pybind11;

std::unique_ptr<SutraGamora> gamora_init(CarmaContext &context, int32_t device,
                                          ArrayFStyle<char> &type, int32_t nactus, int32_t nmodes,
                                          int32_t niter, ArrayFStyle<float> &IFvalue,
                                          ArrayFStyle<int32_t> &IFrowind, ArrayFStyle<int32_t> &IFcolind,
                                          int32_t IFnz, ArrayFStyle<float> &d_TT, ArrayFStyle<float> &pupil,
                                          int32_t size, int32_t Npts, float scale,
                                          ArrayFStyle<float> &Btt, ArrayFStyle<float> &covmodes) {
  return std::unique_ptr<SutraGamora>(new SutraGamora(
      &context, device, type.mutable_data(), nactus, nmodes, niter,
      IFvalue.mutable_data(), IFrowind.mutable_data(), IFcolind.mutable_data(), IFnz,
      d_TT.mutable_data(), pupil.mutable_data(), size, Npts, scale,
      Btt.mutable_data(), covmodes.mutable_data()));
};

int32_t psf_rec_roket(SutraGamora &sg, ArrayFStyle<float> &err) {
  return sg.psf_rec_roket(err.mutable_data());
};

void declare_gamora(py::module &mod) {
  py::class_<SutraGamora>(mod, "Gamora")
      .def(py::init(&gamora_init), R"pbdoc(
    Initializes Gamora

    Args:
          context: (CarmaContext): context

          device: (int): context active device

          type : (str) : reconstruction method used ("roket" or "Vii")

          nactus : (int32_t) : number of actuators

          nmodes (int32_t) : number of modes

          niter : (int32_t) : number of iterations performed with roket

          IFvalue : (np.ndarray[ndim=1,dtype=float32_t]) : Non zeros values of pzt influence function matrix

          IFrowind : (np.ndarray[ndim=1,dtype=int32_t]) : Row indices of nnz values (csr sparse format)

          IFcolind : (np.ndarray[ndim=1,dtype=int32_t]) : Column indices of nnz values (csr sparse format)

          IFnz: (int): number of non zero element in IF

          TT : (np.ndarray[ndim=1,dtype=float32_t])np.ndarray[ndim=1,dtype=float32_t]) : Tip-tilt influence functions

          spupil : (np.ndarray[ndim=2,dtype=float32_t]) : Small pupil

          size: (int): pupil size

          Npts: (int): number of points in the pupil

          scale : (float) : 2*pi/lambda_target with lambda_target expressed in microns

          Btt : (np.ndarray[ndim=2, dtype=np.float32_t]) : Volts to Btt modes matrix

          covmodes : (np.ndarray[ndim=2, dtype=np.float32_t]) : error covariance matrix expressed in a modal basis
           )pbdoc",
           py::arg("context"), py::arg("device"), py::arg("type"),
           py::arg("nactus"), py::arg("nmodes"), py::arg("niter"),
           py::arg("IFvalue"), py::arg("IFrowind"), py::arg("IFcolind"),
           py::arg("IFnz"), py::arg("TT"), py::arg("spupil"), py::arg("size"),
           py::arg("Npts"), py::arg("scale"), py::arg("Btt"),
           py::arg("covmodes"))

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //

      .def_property_readonly(
          "device", [](SutraGamora &sg) { return sg.device; }, "GPU device")

      .def_property_readonly(
          "nactus", [](SutraGamora &sg) { return sg.nactus; },
          "Number of actuators")

      .def_property_readonly(
          "niter", [](SutraGamora &sg) { return sg.niter; },
          "number of iterations")

      .def_property_readonly(
          "nmodes", [](SutraGamora &sg) { return sg.nmodes; },
          "number of modes")

      .def_property_readonly(
          "d_err", [](SutraGamora &sg) { return sg.d_err; }, "Error buffer")

      .def_property_readonly(
          "d_amplipup", [](SutraGamora &sg) { return sg.d_amplipup; },
          "Complex amplitude in the pupil")

      .def_property_readonly(
          "d_psf", [](SutraGamora &sg) { return sg.d_psf; },
          "Reconstructed PSF")

      .def_property_readonly(
          "d_phase", [](SutraGamora &sg) { return sg.d_phase; },
          "Residual phase")

      .def_property_readonly(
          "d_wherephase", [](SutraGamora &sg) { return sg.d_wherephase; },
          "index of valid point")

      .def_property_readonly(
          "d_IF", [](SutraGamora &sg) { return sg.d_IF; }, "sparse IF matrix")

      .def_property_readonly(
          "d_TT", [](SutraGamora &sg) { return sg.d_TT; },
          "tip-tilt IF matrix")

      .def_property_readonly(
          "scale", [](SutraGamora &sg) { return sg.scale; }, "Scale factor")

      .def_property_readonly(
          "size", [](SutraGamora &sg) { return sg.size; },
          "Pupil support size")

      .def_property_readonly(
          "Npts", [](SutraGamora &sg) { return sg.Npts; },
          "number of points in the pupil")

      .def_property_readonly(
          "d_term1", [](SutraGamora &sg) { return sg.d_term1; },
          "Buffer for Vii computation")

      .def_property_readonly(
          "d_term2", [](SutraGamora &sg) { return sg.d_term2; },
          "Buffer for Vii computation")

      .def_property_readonly(
          "d_otftel", [](SutraGamora &sg) { return sg.d_otftel; },
          "OTF of the telescope")

      .def_property_readonly(
          "d_otfVii", [](SutraGamora &sg) { return sg.d_otfVii; },
          "OTF reconstructed from Vii")

      .def_property_readonly(
          "d_mask", [](SutraGamora &sg) { return sg.d_mask; }, "Mask")

      .def_property_readonly(
          "d_eigenvals", [](SutraGamora &sg) { return sg.d_eigenvals; },
          "Eigenvalues of Vii diago")

      .def_property_readonly(
          "d_Btt", [](SutraGamora &sg) { return sg.d_Btt; }, "Btt modal basis")

      .def_property_readonly(
          "d_covmodes", [](SutraGamora &sg) { return sg.d_covmodes; },
          "error covariance marix on the modes")

      .def_property_readonly(
          "d_newmodek", [](SutraGamora &sg) { return sg.d_newmodek; },
          "Mode k from Vii")

      .def_property_readonly(
          "d_Dphi", [](SutraGamora &sg) { return sg.d_Dphi; },
          "Structure function")

      .def_property_readonly(
          "d_pupfft", [](SutraGamora &sg) { return sg.d_pupfft; },
          "FFT of the pupil")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝
      .def("psf_rec_roket", &psf_rec_roket, R"pbdoc(
    Reconstruct the PSF from ROKET error buffer

    Args:
        err: (np.array[ndim=2,dtype=np.float32]): ROKET error buffer

    )pbdoc",
           py::arg("err"))

      .def("psf_rec_Vii", &SutraGamora::psf_rec_Vii, "Vii PSF reconstruction");
};
