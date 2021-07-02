from conans import ConanFile, CMake, tools

import subprocess as sp
from packaging import version

def cuda_version():
    cmd_output = sp.run(['nvcc', '--version'], stdout=sp.PIPE).stdout.decode('utf-8')
    return version.parse(cmd_output.split('release ')[1].split(',')[0])


class CompassConan(ConanFile):
    name = 'compass'
    version = '5.1.0'
    author = 'COMPASS Team <https://github.com/ANR-COMPASS>'
    url = 'https://anr-compass.github.io/compass/'
    description = 'End-to-end AO simulation tool using GPU acceleration'
    topics = ('Adaptive Optics', 'Simulation')
    settings = 'os', 'compiler', 'build_type', 'arch'
    requires = [
        'wyrm/0.4@cosmic/stable',
        'emu/0.1@cosmic/stable',
        ]
    generators = ['cmake', 'cmake_find_package']
    options = {
        'shared'        : [True, False],
        'fPIC'          : [True, False],
        'python_build'  : [True, False],
        'half'          : [True, False],
        'python_version': 'ANY'
    }
    default_options = {
        'shared'      : True,
        'fPIC'        : True,
        'python_build': True,
        'half'        : True,
        'wyrm:cuda'   : True,
        'wyrm:half'   : True,
    }

    def requirements(self):
        if cuda_version() < version.parse('11.0'):
            self.requires('cub/1.8.0@cosmic/stable')

    def _configure(self):
        cmake = CMake(self)
        cmake.definitions['do_half']                 = self.options.half
        cmake.definitions['build_python_module']     = self.options.python_build
        cmake.definitions['PYBIND11_PYTHON_VERSION'] = self.options.python_version
        cmake.definitions['CMAKE_CUDA_ARCHITECTURES'] = self.options['emu'].cuda_sm

        cmake.configure(source_folder='.')

        return cmake

    def build(self):
        self._configure().build()

    def package(self):
        self._configure().install()
