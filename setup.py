import os
import platform
import re
import subprocess
import sys
from distutils.version import LooseVersion

from setuptools import Extension, setup, find_packages
from setuptools.command.build_ext import build_ext

import multiprocessing

from shesha import __version__ as shesha_version

# import shesha
# shesha_version = shesha.shesha.__version__


class CMakeExtension(Extension):

    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuildExt(build_ext):

    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError(
                    "CMake must be installed to build the following extensions: " +
                    ", ".join(e.name for e in self.extensions))

        if platform.system() == "Windows":
            cmake_version = LooseVersion(
                    re.search(r'version\s*([\d.]+)', out.decode()).group(1))
            if cmake_version < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):

        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = [
                '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                '-DPYTHON_EXECUTABLE=' + sys.executable,
                # '-DCMAKE_INSTALL_PREFIX:PATH=.'
        ]

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        if platform.system() == "Windows":
            cmake_args += [
                    '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir)
            ]
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['--', '-j%d' % multiprocessing.cpu_count()]
            # build_args += ['--', 'VERBOSE=1']

        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        if "CUDA_ROOT" in os.environ:
            if os.path.isfile('{}/bin/gcc'.format(os.environ["CUDA_ROOT"])):
                cmake_args += [
                        '-DCMAKE_C_COMPILER={}/bin/gcc'.format(os.environ["CUDA_ROOT"])
                ]
            if os.path.isfile('{}/bin/g++'.format(os.environ["CUDA_ROOT"])):
                cmake_args += [
                        '-DCMAKE_CXX_COMPILER={}/bin/g++'.format(os.environ["CUDA_ROOT"])
                ]

        cmake_args += ['-DVERSION_INFO={}'.format(self.distribution.get_version())]

        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp)
        subprocess.check_call(['cmake', '--build', '.'] + build_args,
                              cwd=self.build_temp)
        # subprocess.check_call(['cmake', '--build', '.', "--target", "install"] + build_args,
        #                       cwd=self.build_temp)


setup(
        name='compass-sim',
        version=shesha_version,
        # author=['Arnaud Sevin'],
        # author_email=['arnaud.sevin@obspm.fr'],
        # description='',
        # long_description='',
        # package_data={'compass-sim': ['lib/pkgconfig/carma.pc', 'lib/pkgconfig/sutra.pc']},
        ext_modules=[CMakeExtension('compass-sim')],
        cmdclass={'build_ext': CMakeBuildExt},
        zip_safe=False,
)

# setup(
#         name='shesha',
#         version='3.0.0',
#         # author=['Arnaud Sevin'],
#         # author_email=['arnaud.sevin@obspm.fr'],
#         # description='',
#         # long_description='',
#         packages=[
#                 'data', 'data.par', 'data.par.par4bench', 'shesha', 'shesha.ao',
#                 'shesha.config', 'shesha.init', 'shesha.scripts', 'shesha.sim',
#                 'shesha.supervisor', 'shesha.util', 'shesha.widgets'
#         ],
#         # packages=find_packages("shesha"),
#         package_dir={
#                 'data': 'shesha/data',
#                 # 'data.par': 'shesha/data/par',
#                 # 'data.par.par4bench': 'shesha/data/par/par4bench',
#                 'shesha': 'shesha/shesha',
#                 # 'shesha.ao': 'shesha/shesha/ao',
#                 # 'shesha.config': 'shesha/shesha/config',
#                 # 'shesha.init': 'shesha/shesha/init',
#                 # 'shesha.scripts': 'shesha/shesha/scripts',
#                 # 'shesha.sim': 'shesha/shesha/sim',
#                 # 'shesha.supervisor': 'shesha/shesha/supervisor',
#                 # 'shesha.sutra_pybind': 'shesha/shesha/sutra_pybind',
#                 # 'shesha.util': 'shesha/shesha/util',
#                 # 'shesha.widgets': 'shesha/shesha/widgets',
#         },
#         package_data={'data': ['layouts/SCAO_PYR.area', 'layouts/SCAO_SH.area']},
#         include_package_data=True,
#         # install_requires=['compass_sim'],
#         zip_safe=False, )
