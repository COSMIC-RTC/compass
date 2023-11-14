from conans import ConanFile, tools, CMake

import os, re, subprocess as sp
from packaging import version

def cuda_version():
    cmd_output = sp.run(['nvcc', '--version'], stdout=sp.PIPE).stdout.decode('utf-8')
    return version.parse(cmd_output.split('release ')[1].split(',')[0])


class CompassConan(ConanFile):
    name = 'compass'
    author = 'COMPASS Team <https://github.com/ANR-COMPASS>'
    url = 'https://anr-compass.github.io/compass/'
    description = 'End-to-end AO simulation tool using GPU acceleration'
    topics = 'Adaptive Optics', 'Simulation'
    settings = 'os', 'compiler', 'build_type', 'arch'
    generators = ['cmake', 'cmake_find_package']

    python_requires = 'cuda_arch/0.1@cosmic/stable'

    options = {
        'shared'        : [True, False],
        'fPIC'          : [True, False],
        'libs'          : [True, False],
        'python'        : [True, False],
        'half'          : [True, False],
        'cuda_sm'       : 'ANY',
        'python_version': 'ANY'
    }
    default_options = {
        'shared'   : True,
        'fPIC'     : True,
        'libs'     : True,
        'python'   : False,
        'python_version': None,
        'half'     : False,
        'wyrm:cuda': True,
        'cuda_sm'  : 'Auto'
    }

    def _cuda_compute_capabilities(self):
        return self.python_requires["cuda_arch"].module.compute_capabilities()

    def _half_support(self):
        archs = str(self.options.cuda_sm)

        # Check if they are all above or equal 60.
        return min(map(int, archs.split(';'))) >= 60

    def set_version(self):
        content = tools.load(os.path.join(self.recipe_folder, 'CMakeLists.txt'))
        version = re.search(r'set\(VERSION_INFO (\d+\.\d+\.\d+)[^\)]*\)', content).group(1)
        self.version = version.strip()

    def configure(self):
        # If `Auto`, replace by local compute capabilities.
        if str(self.options.cuda_sm) == 'Auto':
            self.options.cuda_sm = self._cuda_compute_capabilities()

        if self.options.half:
            self.options.half = self._half_support()

        if self.options.python_version:
            self.options['wyrm'].half = self.options.half


    def requirements(self):
        if cuda_version() < version.parse('11.0'):
            self.requires('cub/1.8.0@cosmic/stable')
        if self.options.python_version:
            self.requires('pybind11/2.10.4')
            self.requires('wyrm/0.4@cosmic/stable')
        else:
            self.options.remove('python_version')

    def _configure(self):
        cmake = CMake(self)
        cmake.definitions['do_half'] = self.options.half
        cmake.definitions['libs_build']     = self.options.libs
        if self.options.python_version:
            cmake.definitions['PYBIND11_PYTHON_VERSION'] = self.options.python_version
            self.options.python = True
        cmake.definitions['python_build']     = self.options.python

        cmake.definitions['CMAKE_CUDA_ARCHITECTURES'] = self.options.cuda_sm

        cmake.configure(source_folder='.')

        return cmake

    def build(self):
        self._configure().build()

    def package(self):
        self._configure().install()
