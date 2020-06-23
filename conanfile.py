from conans import ConanFile, CMake
import subprocess as sp
from packaging import version


def cuda_version():
    cmd_output = sp.run(['nvcc', '--version'], stdout=sp.PIPE).stdout.decode('utf-8')
    return version.parse(cmd_output.split('release ')[1].split(',')[0])


class CompassConan(ConanFile):
    settings = 'os', 'compiler', 'build_type', 'arch'
    requires = ['wyrm/0.1@cosmic/stable']
    generators = 'cmake'
    default_options = {'emu:cuda': True, 'wyrm:cuda': True, 'wyrm:half': True}

    def requirements(self):
        if cuda_version() < version.parse('11.0'):
            self.requires('cub/1.8.0@cosmic/stable')
