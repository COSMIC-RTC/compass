import os
import platform
import re
import subprocess
import sys
from distutils.version import LooseVersion

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext

import multiprocessing


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
        if ext.name != 'naga':
            return

        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))

        cmake_args = [
                '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                '-DPYTHON_EXECUTABLE=' + sys.executable
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

        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        if "CUDA_ROOT" in os.environ:
            if os.path.isfile('{}/bin/gcc'.format(os.environ["CUDA_ROOT"])):
                cmake_args += [
                        '-DCMAKE_C_COMPILER={}/bin/gcc'.format(os.environ["CUDA_ROOT"])
                ]
            if os.path.isfile('{}/bin/g++'.format(os.environ["CUDA_ROOT"])):
                cmake_args += [
                        '-DCMAKE_CXX_COMPILER={}/bin/g++'.format(
                                os.environ["CUDA_ROOT"])
                ]

        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp)
        subprocess.check_call(['cmake', '--build', '.'] + build_args,
                              cwd=self.build_temp)


ext_modules = [
        CMakeExtension('libcarma'),
        CMakeExtension('naga'),
        CMakeExtension('libsutra'),
        CMakeExtension('shesha'),
]

setup(
        name='compass',
        version='1.0.0',
        # author=['Arnaud Sevin'],
        # author_email=['arnaud.sevin@obspm.fr'],
        # description='',
        # long_description='',
        ext_modules=ext_modules,
        cmdclass={'build_ext': CMakeBuildExt},
        zip_safe=False, )

setup(
        name='shesha',
        version='1.0.0',
        # author=['Arnaud Sevin'],
        # author_email=['arnaud.sevin@obspm.fr'],
        # description='',
        # long_description='',
        packages=[
                'shesha', 'shesha.ao', 'shesha.config', 'shesha.init', 'shesha.scripts',
                'shesha.sim', 'shesha.supervisor', 'shesha.sutra_pybind', 'shesha.util',
                'shesha.widgets'
        ],
        package_dir={
                'shesha': 'shesha/shesha',
                'shesha.ao': 'shesha/shesha/ao',
                'shesha.config': 'shesha/shesha/config',
                'shesha.init': 'shesha/shesha/init',
                'shesha.scripts': 'shesha/shesha/scripts',
                'shesha.sim': 'shesha/shesha/sim',
                'shesha.supervisor': 'shesha/shesha/supervisor',
                'shesha.sutra_pybind': 'shesha/shesha/sutra_pybind',
                'shesha.util': 'shesha/shesha/util',
                'shesha.widgets': 'shesha/shesha/widgets',
        },
        install_requires=['compass'],
        zip_safe=False, )
