// -----------------------------------------------------------------------------
//  This file is part of COMPASS <https://anr-compass.github.io/compass/>
//
//  Copyright (C) 2011-2022 COMPASS Team <https://github.com/ANR-COMPASS>
//  All rights reserved.

// -----------------------------------------------------------------------------

//! \file      dms.cpp
//! \ingroup   libsutra
//! \brief     this file provides pybind wrapper for SutraDms
//! \author    COMPASS Team <https://github.com/ANR-COMPASS>
//! \version   5.3.0
//! \date      2022/01/24

#include <sutra_dm.h>

#include <wyrm>

namespace py = pybind11;

std::unique_ptr<SutraDms> dms_init() {
  return std::unique_ptr<SutraDms>(new SutraDms());
};

//  ██████╗ ███╗   ███╗███████╗
//  ██╔══██╗████╗ ████║██╔════╝
//  ██║  ██║██╔████╔██║███████╗
//  ██║  ██║██║╚██╔╝██║╚════██║
//  ██████╔╝██║ ╚═╝ ██║███████║
//  ╚═════╝ ╚═╝     ╚═╝╚══════╝
//

void declare_dms(py::module &mod) {
  py::class_<SutraDms>(mod, "Dms")
      .def(py::init(wy::colCast(dms_init)), R"pbdoc(
            Create a void DMS object
        )pbdoc")

      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //

      .def_property_readonly(
          "d_dms",
          [](SutraDms &sdms) -> vector<SutraDm *> & { return sdms.d_dms; },
          "Vector of SutraDm")

      .def_property_readonly(
          "ndm", [](SutraDms &sdms) { return sdms.d_dms.size(); },
          "Number of SutraDm in SutraDms")
      .def_property_readonly(
          "nact_total", [](SutraDms &sdms) { return sdms.nact_total(); },
          "Total number of actuators in SutraDms")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝

      .def("add_dm",
           (int (SutraDms::*)(CarmaContext *, const char *, float, long, long,
                              long, long, long, float, long, float, float,
                              float, float, int)) &
               SutraDms::add_dm,
           R"pbdoc(
    Add a SutraDm in the SutraDms vector

    Args:
      context: (CarmaContext) : current carma context

      type: (str): DM type ("pzt", "kl", or "tt")

      alt: (float): Conjugaison altitude in meters

      dim: (long): Support dimension

      nactus: (long): Number of actuators

      influsize: (long): Influenction function support size

      ninflupos: (long): Size of _influpos array

      n_npoints: (long): Size of _ninflu array

      push4imat: (float): Voltage to apply for imat computation

      nord: (long): Number of radial order for kl dm (0 if not kl)

      dx: (float): X axis misregistration [pixels]

      dy: (float): Y axis misregistration [pixels]

      theta: (float): Rotation angle misregistration [radians]

      G: (float): Magnification factor

      device: (int): Device index
      )pbdoc",
           py::arg("context"), py::arg("type"), py::arg("alt"), py::arg("dim"),
           py::arg("nactus"), py::arg("influsize"), py::arg("ninflupos"),
           py::arg("n_npoints"), py::arg("push4imat"), py::arg("nord"),
           py::arg("dx"), py::arg("dy"), py::arg("thetaML"), py::arg("G"),
           py::arg("device"))

      .def("add_dm",
           (int (SutraDms::*)(CarmaContext *, const char *, float, long, long,
                              long, long, long, float, long, int)) &
               SutraDms::add_dm,
           R"pbdoc(
    Add a SutraDm in the SutraDms vector

    Args:
      context: (CarmaContext) : current carma context

      type: (str): DM type ("pzt", "kl", or "tt")

      alt: (float): Conjugaison altitude in meters

      dim: (long): Support dimension

      nactus: (long): Number of actuators

      influsize: (long): Influenction function support size

      ninflupos: (long): Size of _influpos array

      n_npoints: (long): Size of _ninflu array

      push4imat: (float): Voltage to apply for imat computation

      nord: (long): Number of radial order for kl dm (0 if not kl)

      device: (int): Device index
        )pbdoc",
           py::arg("context"), py::arg("type"), py::arg("alt"), py::arg("dim"),
           py::arg("nactus"), py::arg("influsize"), py::arg("ninflupos"),
           py::arg("n_npoints"), py::arg("push4imat"), py::arg("nord"),
           py::arg("device"))

      .def("insert_dm", wy::colCast(&SutraDms::insert_dm),
           R"pbdoc(
    Add a SutraDm in the SutraDms vector at the specified index

    Args:
      context: (CarmaContext) : current carma context

      type: (str): DM type ("pzt", "kl", or "tt")

      alt: (float): Conjugaison altitude in meters

      dim: (long): Support dimension

      nactus: (long): Number of actuators

      influsize: (long): Influenction function support size

      ninflupos: (long): Size of _influpos array

      n_npoints: (long): Size of _ninflu array

      push4imat: (float): Voltage to apply for imat computation

      nord: (long): Number of radial order for kl dm (0 if not kl)

      dx: (float): X axis misregistration [pixels]

      dy: (float): Y axis misregistration [pixels]

      theta: (float): Rotation angle misregistration [radians]

      G: (float): Magnification factor

      device: (int): Device index

      idx: (int) : DM index in the vector dms
      )pbdoc",
           py::arg("context"), py::arg("type"), py::arg("alt"), py::arg("dim"),
           py::arg("nactus"), py::arg("influsize"), py::arg("ninflupos"),
           py::arg("n_npoints"), py::arg("push4imat"), py::arg("nord"),
           py::arg("dx"), py::arg("dy"), py::arg("theta"), py::arg("G"),
           py::arg("device"), py::arg("idx"))

      .def("remove_dm", wy::colCast(&SutraDms::remove_dm),
           R"pbdoc(
    Remove and delete the selected DM from SutraDms

    Args:
      idx: (int): index of DM
      )pbdoc",
           py::arg("idx"))

      .def("__str__",
           [](SutraDms &sdms) {
             std::cout << sdms.d_dms.size() << " DMs created" << std::endl;
             std::cout << "Total number of actuators : " << sdms.nact_total()
                       << std::endl;
             std::cout << "DM # | Type  |   Alt   | Nact | Dim" << std::endl;
             vector<SutraDm *>::iterator it = sdms.d_dms.begin();
             int i = 0;
             while (it != sdms.d_dms.end()) {
               std::cout << i << " | " << (*it)->type << " | "
                         << (*it)->altitude << " | " << (*it)->nactus << " | "
                         << (*it)->dim << std::endl;
               i++;
               it++;
             }
             return "";
           })

      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
      //

      .def(
          "set_full_com",
          [](SutraDms &sdms,
             py::array_t<float, py::array::f_style | py::array::forcecast> data,
             bool shape_dm = true) {
            if (data.size() != sdms.nact_total())
              DEBUG_TRACE("Command vector has wrong dimension");
            else {
              int com_idx = 0;
              float *com = data.mutable_data();
              for (vector<SutraDm *>::iterator it = sdms.d_dms.begin();
                   it != sdms.d_dms.end(); it++) {
                (*it)->d_com->host2device(&com[com_idx]);
                com_idx += (*it)->nactus;
                if (shape_dm) (*it)->comp_shape();
              }
            }
          },
          R"pbdoc(
    Set the command vector of all DM in SutraDms, and computes the DMs shapes

    Args:
        com: (np.array(ndim=1, dtype=np.float32)): Concatened command vectors

        shape_dm: (bool): (optionnal, default=True)  Computes the DM shape
      )pbdoc",
          py::arg("com"), py::arg("shape_dm") = true)

      ;
};

//  ██████╗ ███╗   ███╗
//  ██╔══██╗████╗ ████║
//  ██║  ██║██╔████╔██║
//  ██║  ██║██║╚██╔╝██║
//  ██████╔╝██║ ╚═╝ ██║
//  ╚═════╝ ╚═╝     ╚═╝
//

void declare_dm(py::module &mod) {
  py::class_<SutraDm>(mod, "Dm")
      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //
      .def_property_readonly(
          "device", [](SutraDm &sdm) { return sdm.device; }, "Device index")

      .def_property_readonly(
          "type", [](SutraDm &sdm) { return sdm.type; }, "DM type")

      .def_property_readonly(
          "altitude", [](SutraDm &sdm) { return sdm.altitude; },
          "DM conjugaison altitude")

      .def_property_readonly(
          "nactus", [](SutraDm &sdm) { return sdm.nactus; },
          "Number of actuators")

      .def_property_readonly(
          "influsize", [](SutraDm &sdm) { return sdm.influsize; },
          "Influence function support size")

      .def_property_readonly(
          "dim", [](SutraDm &sdm) { return sdm.dim; }, "DM support size")

      .def_property_readonly(
          "push4imat", [](SutraDm &sdm) { return sdm.push4imat; },
          "Voltage to apply for imat computation")

      .def_property_readonly(
          "d_shape", [](SutraDm &sdm) { return sdm.d_shape->d_screen; },
          "DM shape")

      .def_property_readonly(
          "d_com", [](SutraDm &sdm) { return sdm.d_com; },
          "Current commands of the DM")

      .def_property_readonly(
          "d_influ", [](SutraDm &sdm) { return sdm.d_influ; },
          "Cube of influence functions")

      .def_property_readonly(
          "d_istart", [](SutraDm &sdm) { return sdm.d_istart; },
          "TODO: docstring")

      .def_property_readonly(
          "d_npoints", [](SutraDm &sdm) { return sdm.d_npoints; },
          "Number of IF that impact each pixel of the pupil")

      .def_property_readonly(
          "d_influpos", [](SutraDm &sdm) { return sdm.d_influpos; },
          "Influence functions positions in the pupil")

      .def_property_readonly(
          "d_xoff", [](SutraDm &sdm) { return sdm.d_xoff; }, "TODO: docstring")

      .def_property_readonly(
          "dx", [](SutraDm &sdm) { return sdm.dx; }, "X registration in pixels")

      .def_property_readonly(
          "dy", [](SutraDm &sdm) { return sdm.dy; }, "Y registration in pixels")

      .def_property_readonly(
          "thetaML", [](SutraDm &sdm) { return sdm.thetaML; },
          "thetaML registration in radians")

      .def_property_readonly(
          "G", [](SutraDm &sdm) { return sdm.G; },
          "Magnification factor registration in pixels")

      .def_property_readonly(
          "d_yoff", [](SutraDm &sdm) { return sdm.d_yoff; }, "TODO: docstring")

      .def_property_readonly(
          "d_KLbasis", [](SutraDm &sdm) { return sdm.d_KLbasis; },
          "KL to volts matrix")

      .def_property_readonly(
          "d_kl", [](SutraDm &sdm) { return sdm.d_kl; }, "SutraKL DM")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝
      .def("pzt_loadarrays", wy::colCast(&SutraDm::pzt_loadarrays),
           R"pbdoc(
    Load all the arrays computed during the initialization for a pzt DM in a SutraDm object

    Args:
        influ: (np.ndarray[ndim=3,dtype=np.float32_t]) : influence functions cube

        influ2:

        influ3:

        influpos: (np.ndarray[ndim=1,dtype=np.int32_t]) : positions of the IF in the pupil

        influpos2:

        npoints: (np.ndarray[ndim=1,dtype=np.int32_t]) : for each pixel on the DM screen, the number of IF which impact on this pixel

        istart: (np.ndarray[ndim=1,dtype=np.int32_t]) :

        xoff: (np.ndarray[ndim=1,dtype=np.int32_t]) : x-offset for shape computation

        yoff: (np.ndarray[ndim=1,dtype=np.int32_t]) : y-offset or shape computation
        )pbdoc",
           py::arg("influ"),
           /*py::arg("influ2"),py::arg("influ3"),*/ py::arg("influpos"),
           /*py::arg("influpos2"),*/ py::arg("npoints"), py::arg("istart"),
           py::arg("xoff"), py::arg("yoff"))

      .def("tt_loadarrays", wy::colCast(&SutraDm::tt_loadarrays),
           R"pbdoc(
    Load all the arrays computed during the initialization for a tt DM in a SutraDm object

    Args:
        influ: (np.ndarray[ndim=3,dtype=np.float32_t]) : influence functions cube
        )pbdoc",
           py::arg("influ"))
      // .def("tt_loadarrays",[](SutraDm &sdm, py::array_t<float,
      // py::array::f_style | py::array::forcecast> data){
      //   sdm.d_influ->host2device(data.mutable_data());
      // })
      .def("kl_loadarrays", wy::colCast(&SutraDm::kl_loadarrays),
           R"pbdoc(
    Load all the arrays computed during the initialization for a kl DM in a SutraDm object

    Args:
        rabas: (np.ndarray[ndim=1,dtype=np.float32_t]): TODO docstring

        azbas: (np.ndarray[ndim=1,dtype=np.float32_t]):

        ords: (np.ndarray[ndim=1,dtype=np.int32_t]):

        cr: (np.ndarray[ndim=1,dtype=np.float32_t]):

        cp: (np.ndarray[ndim=1,dtype=np.float32_t]):
        )pbdoc",
           py::arg("rabas"), py::arg("azbas"), py::arg("ords"), py::arg("cr"),
           py::arg("cp"))

      .def("reset_shape", &SutraDm::reset_shape, R"pbdoc(
        Reset the DM shape to zeros (flat)
      )pbdoc")

      .def("comp_shape", (int (SutraDm::*)(void)) & SutraDm::comp_shape,
           R"pbdoc(
        Compute the DM shape according to commands d_com
      )pbdoc")

      .def("comp_shape",
           wy::colCast((int (SutraDm::*)(float *)) & SutraDm::comp_shape),
           R"pbdoc(
    Compute the DM shape according to given commands

    Args:
        com: (np.ndarray[ndim=1,dtype=np.float32_t]): commands to apply

      )pbdoc",
           py::arg("com"))

      .def("comp_shape",
           wy::colCast((int (SutraDm::*)(uint16_t *)) & SutraDm::comp_shape),
           R"pbdoc(
    Compute the DM shape according to given commands

    Args:
        com: (np.ndarray[ndim=1,dtype=np.uint16_t]): commands to apply

      )pbdoc",
           py::arg("com"))

      .def("comp_oneactu", wy::colCast(&SutraDm::comp_oneactu), R"pbdoc(
    Push the specified actuator and computes the corresponding DM shape

    Args:
        nactu: (int): Actuator index

        ampli: (float): Volt to apply to this actuator
      )pbdoc",
           py::arg("nactu"), py::arg("ampli"))

      .def("set_registration", wy::colCast(&SutraDm::set_registration), R"pbdoc(
    Set the registration parameters : dx, dy, theta and G

    Args:
        dx: (float): X axis misregistration [pixels]

        dy: (float): Y axis misregistration [pixels]

        theta: (float): Rotation angle misregistration [radians]

        G: (float): Magnification factor
      )pbdoc",
           py::arg("dx"), py::arg("dy"), py::arg("theta"), py::arg("G"))

      .def("compute_KLbasis", wy::colCast(&SutraDm::compute_KLbasis),
           R"pbdoc(
    Computes the KL to volt matrix by double diagonalisation (cf. Gendron thesis)
        - compute the phase covariance matrix on the actuators using Kolmogorov
        - compute the geometric covariance matrix
        - double diagonalisation to obtain KL basis

    Args:
        xpos: (np.ndarray[ndim=1,dtype=np.float32_t]) : x-position of actuators

        ypos: (np.ndarray[ndim=1,dtype=np.float32_t]) : y-position of actuators

        indx_pup: (np.ndarray[ndim=1,dtype=np.int32_t]) : indices of pupil points

        dim: (long) : number of points in the pupil

        norm: (float) : normalization factor

        ampli: (float) : amplitude
      )pbdoc",
           py::arg("xpos"), py::arg("ypos"), py::arg("indx_pup"),
           py::arg("dim"), py::arg("norm"), py::arg("ampli"))

      .def("__str__",
           [](SutraDm &sdm) {
             std::cout << "Type  |   Alt   | Nact | Dim" << std::endl;
             std::cout << sdm.type << " | " << sdm.altitude << " | "
                       << sdm.nactus << " | " << sdm.dim << std::endl;
             return "";
           })

      //  ███████╗███████╗████████╗████████╗███████╗██████╗ ███████╗
      //  ██╔════╝██╔════╝╚══██╔══╝╚══██╔══╝██╔════╝██╔══██╗██╔════╝
      //  ███████╗█████╗     ██║      ██║   █████╗  ██████╔╝███████╗
      //  ╚════██║██╔══╝     ██║      ██║   ██╔══╝  ██╔══██╗╚════██║
      //  ███████║███████╗   ██║      ██║   ███████╗██║  ██║███████║
      //  ╚══════╝╚══════╝   ╚═╝      ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝
      //

      .def(
          "set_com",
          [](SutraDm &sdm,
             py::array_t<float, py::array::f_style | py::array::forcecast> data,
             bool shape_dm = true) {
            if (sdm.d_com->get_nb_elements() == data.size()) {
              sdm.d_com->host2device(data.mutable_data());
              if (shape_dm) sdm.comp_shape();
            } else
              DEBUG_TRACE("Wrong dimension");
          },
          R"pbdoc(
    Set the command vector of a SutraDm, and computes the DM shape

    Args:
        com: (np.array(ndim=1, dtype=np.float32)): Command vector

        shape_dm: (bool): (optionnal, default=True)  Computes the DM shape
      )pbdoc",
          py::arg("com"), py::arg("shape_dm") = true);
};
