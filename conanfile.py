from conans import ConanFile, CMake, tools

import os
import subprocess as sp
from packaging import version


def cuda_version():
    cmd_output = sp.run(['nvcc', '--version'], stdout=sp.PIPE).stdout.decode('utf-8')
    return version.parse(cmd_output.split('release ')[1].split(',')[0])


class CompassConan(ConanFile):
    name = "compass"
    version = "5.0"
    author = "COMPASS Team <https://github.com/ANR-COMPASS>"
    url = "https://anr-compass.github.io/compass/"
    description = "End-to-end AO simulation tool using GPU acceleration"
    topics = ("Adaptive Optics", "Simulation")
    settings = 'os', 'compiler', 'build_type', 'arch'
    requires = ['wyrm/0.1@cosmic/stable']
    generators = 'cmake'
    default_options = {'wyrm:cuda': True, 'wyrm:half': False}
    scm = {
            "type": "git",
            "subfolder": ".",
            "url": "https://gitlab.obspm.fr/compass/compass.git",
            "revision": "rc",
            "username": "conan",
            "password": "yLQJgzs1-x_yYWxjWerA"
    }

    def requirements(self):
        if cuda_version() < version.parse('11.0'):
            self.requires('cub/1.8.0@cosmic/stable')

    # def source(self):
    #     self.run("git clone https://gitlab.obspm.fr/compass/compass.git -b rc .")

    def build(self):
        cmake = CMake(self)
        if "COMPASS_DO_HALF" in os.environ:
            cmake.definitions["do_half"] = os.environ['COMPASS_DO_HALF']
        # cmake.configure(source_folder="compass")
        cmake.configure()
        cmake.build()

        # Explicit way:
        # self.run('cmake %s/hello %s'
        #          % (self.source_folder, cmake.command_line))
        # self.run("cmake --build . %s" % cmake.build_config)

    def package(self):
        self.copy("*.h", dst="include", keep_path=False)
        self.copy("lib*/lib*.so*", dst="lib", keep_path=False)
        self.copy("lib*/*.pc", dst="lib/pkgconfig", keep_path=False)
        self.copy("lib*/python_module/*.so", dst="python", keep_path=False)

    #     self.copy("naga", dst="python", keep_path=True)
    #     self.copy("shesha", dst="python", keep_path=True)

    def package_info(self):
        self.cpp_info.libs = ["carma", "sutra"]
