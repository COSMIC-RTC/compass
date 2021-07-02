from conans import ConanFile, tools
from conan.tools.cmake import CMakeToolchain, CMakeDeps, CMake

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

    requires = ['emu/0.1@cosmic/stable']
    python_requires = 'cuda_arch/0.1@cosmic/stable'

    # generators = 'cmake_paths' #, 'cmake_find_package'
    options = {
        'shared'        : [True, False],
        'fPIC'          : [True, False],
        'python'        : [True, False],
        'half'          : [True, False],
        'python_version': 'ANY'
    }
    default_options = {
        'shared'   : True,
        'fPIC'     : True,
        'python'   : True,
        'half'     : True,
        'wyrm:cuda': True,
    }

    def _cuda_compute_capabilities(self):
        return self.python_requires["cuda_arch"].module.compute_capabilities()

    def _half_support(self):
        # Retreive list of cuda compute capabilities.
        if self.options['emu'].cuda_sm == 'Auto':
            archs = self._cuda_compute_capabilities()
        else:
            archs = str(self.options['emu'].cuda_sm)

        # Check if they are all above or equal 60.
        return min(map(int, archs.split(';'))) >= 60

    def set_version(self):
        content = tools.load(os.path.join(self.recipe_folder, 'CMakeLists.txt'))
        version = re.search(r'set\(VERSION_INFO (\d+\.\d+\.\d+)[^\)]*\)', content).group(1)
        self.version = version.strip()

    def configure(self):
        self.options.half = self._half_support()

        if self.options.half and self.options.python:
            self.options['wyrm'].half = self.options.half

    def requirements(self):
        if cuda_version() < version.parse('11.0'):
            self.requires('cub/1.8.0@cosmic/stable')
        if self.options.python:
            self.requires('wyrm/0.4@cosmic/stable')
        else:
            self.options.remove('python_version')

    def generate(self):
        tc = CMakeToolchain(self)
        tc.variables['do_half'] = self.options.half
        tc.variables['build_python_module']     = self.options.python
        if self.options.python:
            tc.variables['PYBIND11_PYTHON_VERSION'] = self.options.python_version

        tc.variables['CMAKE_CUDA_ARCHITECTURES'] = self.options['emu'].cuda_sm
        tc.generate()

        deps = CMakeDeps(self)
        deps.generate()

    def _configure(self):
        cmake = CMake(self)
        cmake.configure()

        return cmake

    def build(self):
        self._configure().build()

    def package(self):
        self._configure().install()
