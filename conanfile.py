from conan import ConanFile
from conan.tools.cmake import CMake, CMakeToolchain, cmake_layout
from conan.tools.files import load

import os
import re


class CompassConan(ConanFile):
    name = "compass"
    author = "COSMIC Team <https://github.com/COSMIC-RTC/compass>"
    url = "https://anr-compass.github.io/compass/"
    description = "End-to-end AO simulation tool using GPU acceleration"
    topics = "Adaptive Optics", "Simulation"
    settings = "os", "compiler", "build_type", "arch"
    generators = ["CMakeDeps"]

    options = {
        "shared":         [True, False],
        "fPIC":           [True, False],
        "libs":           [True, False],
        "half":           [True, False],
        "cuda_sm":        ["native", "all", "all-major", "60", "61", "62",
                           "70", "72", "75", "80", "86", "89", "90"],
        "python_version": [None, "ANY"]
    }
    options_description = {
        "shared":         "Build shared libraries",
        "fPIC":           "Build with -fPIC",
        "libs":           "Build libraries",
        "half":           "Build with half support",
        "cuda_sm":        ("CUDA compute capabilities, native for local "
                           "compute capabilities: native, all, all-major, ",
                           "60, 61, 62, 70, 72, 75, ...\n",
                           "see https://developer.nvidia.com/cuda-gpus for ",
                           "more information"),
        "python_version": "Python version, None for no python bindings",
    }
    default_options = {
        "shared":         True,
        "fPIC":           True,
        "libs":           True,
        "half":           False,
        "cuda_sm":        "native",
        "python_version": None,
    }

    def set_version(self):
        content = load(self,
                       os.path.join(self.recipe_folder, "CMakeLists.txt"))
        version = re.search(r"set\(VERSION_INFO (\d+\.\d+\.\d+)[^\)]*\)",
                            content
                            ).group(1)
        self.version = version.strip()

    def layout(self):
        cmake_layout(self)

    def requirements(self):
        if self.options.python_version:
            self.requires("pybind11/2.11.1")

    def generate(self):
        tc = CMakeToolchain(self)

        print("CUDA compute capabilities: ", self.options.cuda_sm)
        tc.variables['CMAKE_CUDA_ARCHITECTURES'] = self.options.cuda_sm

        # If cuda_sm is not above or equal to 60, disable it.
        tc.variables["do_half"] = self.options.half
        tc.variables["libs_build"] = self.options.libs

        if self.options.python_version:
            tc.variables["PYBIND11_PYTHON_VERSION"] = self.options.python_version
            tc.variables["python_build"] = True
        else:
            tc.variables["python_build"] = False

        tc.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()
        cmake.test()

    def package(self):
        cmake = CMake(self)
        cmake.install()
