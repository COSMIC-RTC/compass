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

//! \file      gamora.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for sutra_gamora
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   4.3.0
//! \date      2011/01/28
//! \copyright GNU Lesser General Public License

#include <wyrm>

#include <sutra_gamora.h>

namespace py = pybind11;
typedef py::array_t<float, py::array::f_style | py::array::forcecast> F_arrayS;

std::unique_ptr<sutra_gamora> gamora_init(carma_context &context, int device,
                                          char *type, int nactus, int nmodes,
                                          int niter, float *IFvalue,
                                          int *IFrowind, int *IFcolind,
                                          int IFnz, float *d_TT, float *pupil,
                                          int size, int Npts, float scale,
                                          float *Btt, float *covmodes) {
  return std::unique_ptr<sutra_gamora>(new sutra_gamora(
      &context, device, type, nactus, nmodes, niter, IFvalue, IFrowind,
      IFcolind, IFnz, d_TT, pupil, size, Npts, scale, Btt, covmodes));
};

void declare_gamora(py::module &mod) {
  py::class_<sutra_gamora>(mod, "Gamora")
      .def(py::init(wy::colCast(gamora_init)), R"pbdoc(
          Initializes Gamora

          Parameters
          ------------
          context: (carma_context): context
          device: (int): context active device
          type : (str) : reconstruction method used ("roket" or "Vii")
          nactus : (int) : number of actuators
          nmodes (int) : number of modes
          niter : (int) : number of iterations performed with roket
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
          "device", [](sutra_gamora &sg) { return sg.device; }, "GPU device")

      .def_property_readonly(
          "nactus", [](sutra_gamora &sg) { return sg.nactus; },
          "Number of actuators")

      .def_property_readonly(
          "niter", [](sutra_gamora &sg) { return sg.niter; },
          "number of iterations")

      .def_property_readonly(
          "nmodes", [](sutra_gamora &sg) { return sg.nmodes; },
          "number of modes")

      .def_property_readonly(
          "d_err", [](sutra_gamora &sg) { return sg.d_err; }, "Error buffer")

      .def_property_readonly(
          "d_amplipup", [](sutra_gamora &sg) { return sg.d_amplipup; },
          "Complex amplitude in the pupil")

      .def_property_readonly(
          "d_psf", [](sutra_gamora &sg) { return sg.d_psf; },
          "Reconstructed PSF")

      .def_property_readonly(
          "d_phase", [](sutra_gamora &sg) { return sg.d_phase; },
          "Residual phase")

      .def_property_readonly(
          "d_wherephase", [](sutra_gamora &sg) { return sg.d_wherephase; },
          "index of valid point")

      .def_property_readonly(
          "d_IF", [](sutra_gamora &sg) { return sg.d_IF; }, "sparse IF matrix")

      .def_property_readonly(
          "d_TT", [](sutra_gamora &sg) { return sg.d_TT; },
          "tip-tilt IF matrix")

      .def_property_readonly(
          "scale", [](sutra_gamora &sg) { return sg.scale; }, "Scale factor")

      .def_property_readonly(
          "size", [](sutra_gamora &sg) { return sg.size; },
          "Pupil support size")

      .def_property_readonly(
          "Npts", [](sutra_gamora &sg) { return sg.Npts; },
          "number of points in the pupil")

      .def_property_readonly(
          "d_term1", [](sutra_gamora &sg) { return sg.d_term1; },
          "Buffer for Vii computation")

      .def_property_readonly(
          "d_term2", [](sutra_gamora &sg) { return sg.d_term2; },
          "Buffer for Vii computation")

      .def_property_readonly(
          "d_otftel", [](sutra_gamora &sg) { return sg.d_otftel; },
          "OTF of the telescope")

      .def_property_readonly(
          "d_otfVii", [](sutra_gamora &sg) { return sg.d_otfVii; },
          "OTF reconstructed from Vii")

      .def_property_readonly(
          "d_mask", [](sutra_gamora &sg) { return sg.d_mask; }, "Mask")

      .def_property_readonly(
          "d_eigenvals", [](sutra_gamora &sg) { return sg.d_eigenvals; },
          "Eigenvalues of Vii diago")

      .def_property_readonly(
          "d_Btt", [](sutra_gamora &sg) { return sg.d_Btt; }, "Btt modal basis")

      .def_property_readonly(
          "d_covmodes", [](sutra_gamora &sg) { return sg.d_covmodes; },
          "error covariance marix on the modes")

      .def_property_readonly(
          "d_newmodek", [](sutra_gamora &sg) { return sg.d_newmodek; },
          "Mode k from Vii")

      .def_property_readonly(
          "d_Dphi", [](sutra_gamora &sg) { return sg.d_Dphi; },
          "Structure function")

      .def_property_readonly(
          "d_pupfft", [](sutra_gamora &sg) { return sg.d_pupfft; },
          "FFT of the pupil")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝
      .def("psf_rec_roket", wy::colCast(&sutra_gamora::psf_rec_roket), R"pbdoc(
        Reconstruct the PSF from ROKET error buffer

        Parameters
        ------------
        err: (np.array[ndim=2,dtype=np.float32]): ROKET error buffer

    )pbdoc",
           py::arg("err"))

      .def("psf_rec_Vii", &sutra_gamora::psf_rec_Vii, "Vii PSF reconstruction");
};
