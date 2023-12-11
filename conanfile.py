from conan import ConanFile
from conan.tools.cmake import CMake, CMakeToolchain, cmake_layout
from conan.tools.files import save, load, get
from conan.tools.layout import basic_layout

import os, re

class CompassConan(ConanFile):
    name = "compass"
    author = "COMPASS Team <https://github.com/ANR-COMPASS>"
    url = "https://anr-compass.github.io/compass/"
    description = "End-to-end AO simulation tool using GPU acceleration"
    topics = "Adaptive Optics", "Simulation"
    settings = "os", "compiler", "build_type", "arch"
    generators = ["CMakeDeps"]

    python_requires = "conan_cuda/1.0.0"

    options = {
        "shared":         [True, False],
        "fPIC":           [True, False],
        "libs":           [True, False],
        "half":           [True, False],
        "cuda_sm":        [None, "ANY"],
        "python_version": [None, "ANY"]
    }
    options_description = {
        "shared":         "Build shared libraries",
        "fPIC":           "Build with -fPIC",
        "libs":           "Build libraries",
        "half":           "Build with half support",
        "cuda_sm":        "CUDA compute capabilities, Auto for local compute capabilities",
        "python_version": "Python version, None for no python bindings"
    }
    default_options = {
        "shared":         True,
        "fPIC":           True,
        "libs":           True,
        "half":           False,
    }

    def set_version(self):
        content = load(self, os.path.join(self.recipe_folder, "CMakeLists.txt"))
        version = re.search(r"set\(VERSION_INFO (\d+\.\d+\.\d+)[^\)]*\)",
                            content
                            ).group(1)
        self.version = version.strip()

    def layout(self):
        cmake_layout(self)

    def requirements(self):
        cuda_major = self.python_requires["conan_cuda"].module.properties().major
        if int(cuda_major) < 11:
            self.requires("cub/1.8.0@cosmic/stable")

        if self.options.python_version:
            self.requires("pybind11/2.11.1")

    def generate(self):
        tc = CMakeToolchain(self)

        # cuda helper
        cuda = self.python_requires['conan_cuda'].module

        # If cuda_sm is not set, detect local compute capabilities.
        sm = self.options.cuda_sm if self.options.cuda_sm else cuda.compute_capabilities()

        print("CUDA compute capabilities: ", sm)
        tc.variables['CMAKE_CUDA_ARCHITECTURES'] = sm

        # Check if they are all above or equal 60. Disable half if not.
        tc.variables["do_half"] = self.options.half and min(map(int, str(sm).split(";"))) >= 60

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
