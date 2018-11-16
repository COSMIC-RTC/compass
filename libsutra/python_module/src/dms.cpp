#include <wyrm>

#include <sutra_dm.h>

namespace py = pybind11;

std::unique_ptr<sutra_dms> dms_init() {
  return std::unique_ptr<sutra_dms>(new sutra_dms());
};

//  ██████╗ ███╗   ███╗███████╗
//  ██╔══██╗████╗ ████║██╔════╝
//  ██║  ██║██╔████╔██║███████╗
//  ██║  ██║██║╚██╔╝██║╚════██║
//  ██████╔╝██║ ╚═╝ ██║███████║
//  ╚═════╝ ╚═╝     ╚═╝╚══════╝
//

void declare_dms(py::module &mod) {
  py::class_<sutra_dms>(mod, "Dms")
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
          [](sutra_dms &sdms) -> vector<sutra_dm *> & { return sdms.d_dms; },
          "Vector of sutra_dm")

      .def_property_readonly("ndm",
                             [](sutra_dms &sdms) { return sdms.d_dms.size(); },
                             "Number of sutra_dm in sutra_dms")
      .def_property_readonly("nact_total",
                             [](sutra_dms &sdms) { return sdms.nact_total(); },
                             "Total number of actuators in sutra_dms")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝

      .def("add_dm", wy::colCast(&sutra_dms::add_dm),
           R"pbdoc(
        Add a sutra_dm in the sutra_dms vector

        Parameters
        ------------
        context: (carma_context) : current carma context
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

      .def("remove_dm", wy::colCast(&sutra_dms::remove_dm),
           R"pbdoc(
        Remove and delete the selected DM from sutra_dms
        Parameters
        ------------
        idx: (int): index of DM
        )pbdoc",
           py::arg("idx"))

      .def("__str__",
           [](sutra_dms &sdms) {
             std::cout << sdms.d_dms.size() << " DMs created" << std::endl;
             std::cout << "Total number of actuators : " << sdms.nact_total()
                       << std::endl;
             std::cout << "DM # | Type  |   Alt   | Nact | Dim" << std::endl;
             vector<sutra_dm *>::iterator it = sdms.d_dms.begin();
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
          [](sutra_dms &sdms,
             py::array_t<float, py::array::f_style | py::array::forcecast> data,
             bool shape_dm = true) {
            if (data.size() != sdms.nact_total())
              DEBUG_TRACE("Command vector has wrong dimension");
            else {
              int com_idx = 0;
              float *com = data.mutable_data();
              for (vector<sutra_dm *>::iterator it = sdms.d_dms.begin();
                   it != sdms.d_dms.end(); it++) {
                (*it)->d_com->host2device(&com[com_idx]);
                com_idx += (*it)->nactus;
                if (shape_dm) (*it)->comp_shape();
              }
            }
          },
          R"pbdoc(
        Set the command vector of all DM in sutra_dms, and computes the DMs shapes

        Parameters
        ------------
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
  py::class_<sutra_dm>(mod, "Dm")
      //  ██████╗ ██████╗  ██████╗ ██████╗ ███████╗██████╗ ████████╗██╗   ██╗
      //  ██╔══██╗██╔══██╗██╔═══██╗██╔══██╗██╔════╝██╔══██╗╚══██╔══╝╚██╗ ██╔╝
      //  ██████╔╝██████╔╝██║   ██║██████╔╝█████╗  ██████╔╝   ██║    ╚████╔╝
      //  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ ██╔══╝  ██╔══██╗   ██║     ╚██╔╝
      //  ██║     ██║  ██║╚██████╔╝██║     ███████╗██║  ██║   ██║      ██║
      //  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═╝   ╚═╝      ╚═╝
      //
      .def_property_readonly("device", [](sutra_dm &sdm) { return sdm.device; },
                             "Device index")

      .def_property_readonly("type", [](sutra_dm &sdm) { return sdm.type; },
                             "DM type")

      .def_property_readonly("altitude",
                             [](sutra_dm &sdm) { return sdm.altitude; },
                             "DM conjugaison altitude")

      .def_property_readonly("nactus", [](sutra_dm &sdm) { return sdm.nactus; },
                             "Number of actuators")

      .def_property_readonly("influsize",
                             [](sutra_dm &sdm) { return sdm.influsize; },
                             "Influence function support size")

      .def_property_readonly("dim", [](sutra_dm &sdm) { return sdm.dim; },
                             "DM support size")

      .def_property_readonly("push4imat",
                             [](sutra_dm &sdm) { return sdm.push4imat; },
                             "Voltage to apply for imat computation")

      .def_property_readonly(
          "d_shape", [](sutra_dm &sdm) { return sdm.d_shape->d_screen; },
          "DM shape")

      .def_property_readonly("d_com", [](sutra_dm &sdm) { return sdm.d_com; },
                             "Current commands of the DM")

      .def_property_readonly("d_influ",
                             [](sutra_dm &sdm) { return sdm.d_influ; },
                             "Cube of influence functions")

      .def_property_readonly("d_istart",
                             [](sutra_dm &sdm) { return sdm.d_istart; },
                             "TODO: docstring")

      .def_property_readonly("d_npoints",
                             [](sutra_dm &sdm) { return sdm.d_npoints; },
                             "Number of IF that impact each pixel of the pupil")

      .def_property_readonly("d_influpos",
                             [](sutra_dm &sdm) { return sdm.d_influpos; },
                             "Influence functions positions in the pupil")

      .def_property_readonly("d_xoff", [](sutra_dm &sdm) { return sdm.d_xoff; },
                             "TODO: docstring")

      .def_property_readonly("d_yoff", [](sutra_dm &sdm) { return sdm.d_yoff; },
                             "TODO: docstring")

      .def_property_readonly("d_KLbasis",
                             [](sutra_dm &sdm) { return sdm.d_KLbasis; },
                             "KL to volts matrix")

      .def_property_readonly("d_kl", [](sutra_dm &sdm) { return sdm.d_kl; },
                             "sutra_kl DM")

      //  ███╗   ███╗███████╗████████╗██╗  ██╗ ██████╗ ██████╗ ███████╗
      //  ████╗ ████║██╔════╝╚══██╔══╝██║  ██║██╔═══██╗██╔══██╗██╔════╝
      //  ██╔████╔██║█████╗     ██║   ███████║██║   ██║██║  ██║███████╗
      //  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██║  ██║╚════██║
      //  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝██████╔╝███████║
      //  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝
      .def("pzt_loadarrays", wy::colCast(&sutra_dm::pzt_loadarrays),
           R"pbdoc(
        Load all the arrays computed during the initialization
        for a pzt DM in a sutra_dm object

        Parameters
        ------------
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

      .def("tt_loadarrays", wy::colCast(&sutra_dm::tt_loadarrays),
           R"pbdoc(
        Load all the arrays computed during the initialization
        for a tt DM in a sutra_dm object

        Parameters
        ------------
        influ: (np.ndarray[ndim=3,dtype=np.float32_t]) : influence functions cube
        )pbdoc",
           py::arg("influ"))
      // .def("tt_loadarrays",[](sutra_dm &sdm, py::array_t<float,
      // py::array::f_style | py::array::forcecast> data){
      //   sdm.d_influ->host2device(data.mutable_data());
      // })
      .def("kl_loadarrays", wy::colCast(&sutra_dm::kl_loadarrays),
           R"pbdoc(
        Load all the arrays computed during the initialization
        for a kl DM in a sutra_dm object

        Parameters
        ------------
        rabas: (np.ndarray[ndim=1,dtype=np.float32_t]): TODO docstring
        azbas: (np.ndarray[ndim=1,dtype=np.float32_t]):
        ords: (np.ndarray[ndim=1,dtype=np.int32_t]):
        cr: (np.ndarray[ndim=1,dtype=np.float32_t]):
        cp: (np.ndarray[ndim=1,dtype=np.float32_t]):
        )pbdoc",
           py::arg("rabas"), py::arg("azbas"), py::arg("ords"), py::arg("cr"),
           py::arg("cp"))

      .def("reset_shape", &sutra_dm::reset_shape, R"pbdoc(
        Reset the DM shape to zeros (flat)
      )pbdoc")

      .def("comp_shape", (int (sutra_dm::*)(void)) & sutra_dm::comp_shape,
           R"pbdoc(
        Compute the DM shape according to commands d_com
      )pbdoc")

      .def("comp_shape",
           wy::colCast((int (sutra_dm::*)(float *)) & sutra_dm::comp_shape),
           R"pbdoc(
        Compute the DM shape according to given commands

        Parameters
        ------------
        com: (np.ndarray[ndim=1,dtype=np.float32_t]): commands to apply

      )pbdoc",
           py::arg("com"))

      .def("comp_oneactu", wy::colCast(&sutra_dm::comp_oneactu), R"pbdoc(
        Push the specified actuator and computes the corresponding DM shape

        Parameters
        ------------
        nactu: (int): Actuator index
        ampli: (float): Volt to apply to this actuator
      )pbdoc",
           py::arg("nactu"), py::arg("ampli"))

      .def("compute_KLbasis", wy::colCast(&sutra_dm::compute_KLbasis),
           R"pbdoc(
        Computes the KL to volt matrix by double diagonalisation (cf. Gendron thesis)
            - compute the phase covariance matrix on the actuators using Kolmogorov
            - compute the geometric covariance matrix
            - double diagonalisation to obtain KL basis

        Parameters
        ------------
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
           [](sutra_dm &sdm) {
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
          [](sutra_dm &sdm,
             py::array_t<float, py::array::f_style | py::array::forcecast> data,
             bool shape_dm = true) {
            if (sdm.d_com->getNbElem() == data.size())
              sdm.d_com->host2device(data.mutable_data());
            else
              DEBUG_TRACE("Wrong dimension");
          },
          R"pbdoc(
        Set the command vector of a sutra_dm, and computes the DM shape

        Parameters
        ------------
        com: (np.array(ndim=1, dtype=np.float32)): Command vector
        shape_dm: (bool): (optionnal, default=True)  Computes the DM shape
      )pbdoc",
          py::arg("com"), py::arg("shape_dm") = true);
};
