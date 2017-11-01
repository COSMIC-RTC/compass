# encoding: utf-8

import sys
import os
from os.path import isfile

import subprocess
from distutils.core import setup

import codecs

# enable to dump annotate html for each pyx source file
# import Cython.Compiler.Options
# Cython.Compiler.Options.annotate = True

# import shutil
from Cython.Build import cythonize

import numpy
from distutils.extension import Extension

# Remove the "-Wstrict-prototypes" compiler option, which isn't valid for C++.
import distutils.sysconfig
cfg_vars = distutils.sysconfig.get_config_vars()
for key, value in list(cfg_vars.items()):
    if type(value) == str:
        cfg_vars[key] = value.replace("-Wstrict-prototypes", "")
# ==================================

print((sys.prefix))
print("SYS.PATH")
print("======================================")
print((sys.path))
print("======================================")

listMod = ["Telescope", "Sensors", "Atmos", "Dms", "Target",
           "Rtc"]  # , "Gamora", "Groot"]
dependencies = {
        "Sensors": ["Telescope"],
        "Target": ["Telescope"]
        # "shesha_roket": listMod[:-2],
}

naga_path = os.environ.get('NAGA_ROOT')
if (naga_path is None):
    raise EnvironmentError("Environment variable 'NAGA_ROOT' must be define")
sys.path.append(naga_path + '/src')

shesha_path = os.environ.get('SHESHA_ROOT')
if (shesha_path is None):
    raise EnvironmentError("Environment variable 'SHESHA_ROOT' must be define")

if not os.path.exists(shesha_path + "/lib"):
    os.makedirs(shesha_path + "/lib")

sys.path.append(shesha_path + "/lib")


def locate_compass():
    """Locate compass library"""

    if 'COMPASS_ROOT' in os.environ:
        root_compass = os.environ['COMPASS_ROOT']

    else:
        raise EnvironmentError("Environment variable 'COMPASS_ROOT' must be define")

    compass_config = {
            'inc_sutra': root_compass + '/libsutra/include.h',
            'inc_carma': root_compass + '/libcarma/include.h',
            'inc_naga': root_compass + '/naga',
            'lib': root_compass
    }

    return compass_config


COMPASS = locate_compass()
# Obtain the numpy include directory.  This logic works across numpy versions.
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

#source = ['shesha']
libraries = ['sutra']
include_dirs = [
        numpy_include, COMPASS['inc_carma'], COMPASS['inc_sutra'], COMPASS['inc_naga']
]  # , shesha_path+"/src"]

library_dirs = [COMPASS['lib'] + "/libsutra"]

if 'CUDA_INC_PATH' in os.environ:
    include_dirs.append(os.environ['CUDA_INC_PATH'])
else:
    raise EnvironmentError("Environment variable 'CUDA_INC_PATH' must be define")

# deprecated
# def which(program):
#     import os
#
#     def is_exe(fpath):
#         return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
#
#     fpath, fname = os.path.split(program)
#     if fpath:
#         if is_exe(program):
#             return program
#     else:
#         for path in os.environ["PATH"].split(os.pathsep):
#             path = path.strip('"')
#             exe_file = os.path.join(path, program)
#             if is_exe(exe_file):
#                 return exe_file
#
#     return None

parFile = None
if not isfile("par.pxi"):
    parFile = codecs.open("par.pxi", mode='w', encoding='utf-8')
else:
    import warnings
    warnings.warn("par.pxi found, it will not be updated", Warning)

USE_MPI = 0

# if which("mpicxx"):
#     try:
#         import mpi4py
#         MPI=locate_MPI()
#         libraries.extend(MPI['clibs'])
#         include_dirs.extend(MPI['cincdirs'])
#         include_dirs.extend(MPI['pincdirs'])
#         library_dirs.extend(MPI['clibdirs'])
#         USE_MPI = 2
#     except ImportError:
#         print("mpi4py not found, MPI disabled")
# else:
#     print("mpicxx not found, MPI disabled")
if parFile:
    parFile.write("DEF USE_MPI=%d # 0/1/2 \n" % USE_MPI)

USE_BRAMA = 0
define_macros = []
if 'BRAMA_ROOT' in os.environ:
    brama_root = os.environ.get('BRAMA_ROOT')
    ace_root = os.environ.get('ACE_ROOT')
    tao_root = os.environ.get('TAO_ROOT')
    dds_root = os.environ.get('DDS_ROOT')
    USE_BRAMA = 1
    include_dirs.extend([brama_root])
    include_dirs.extend([ace_root])
    include_dirs.extend([tao_root])
    include_dirs.extend([tao_root + '/orbsvcs'])
    include_dirs.extend([dds_root])
    library_dirs.extend([ace_root + '/lib'])
    library_dirs.extend([tao_root + '/tao'])
    library_dirs.extend([dds_root + '/lib'])
    library_dirs.extend([brama_root])
    libraries.extend(['BRAMACommon'])
    libraries.extend(['OpenDDS_InfoRepoDiscovery'])
    libraries.extend(['OpenDDS_Dcps'])
    libraries.extend(['TAO_PortableServer'])
    libraries.extend(['TAO_AnyTypeCode'])
    libraries.extend(['TAO'])
    libraries.extend(['ACE'])
    libraries.extend(['dl'])
    libraries.extend(['rt'])
    define_macros = [
            ('USE_BRAMA', None),
            ('_GNU_SOURCE', None),
            ('__ACE_INLINE__', None),
    ]

if parFile:
    parFile.write("DEF USE_BRAMA=%d # 0/1 \n" % USE_BRAMA)
    parFile.close()


def dependencies_module(name):
    print("=======================================")
    print("resolving dependencies for", name)
    print("=======================================")
    try:
        dep = dependencies[name]
        print(("dependencies:", dep))
        if (os.path.exists("src/sutra_bind/" + name + ".cpp")):
            for d in dep:
                if (os.stat("src/sutra_bind/" + d + ".pyx").st_mtime >
                            os.stat("src/sutra_bind/" + name + ".cpp").st_mtime):
                    # cpp file outdated if exists
                    if (os.path.exists("src/sutra_bind/" + name + ".cpp")):
                        os.remove("src/sutra_bind/" + name + ".cpp")
    except:  # KeyError e:
        print("No depencies found")


def compile_module(name):
    print("=======================================")
    print(("creating module ", name))
    print("=======================================")
    ext = Extension(
            shesha_path + "/lib/" + name,
            sources=['src/sutra_bind/' + name + '.pyx'],
            extra_compile_args=[
                    "-Wno-unused-function",
                    "-Wno-unused-label",
                    "-Wno-cpp",
                    "-Wno-deprecated-declarations",
                    "-std=c++11",
                    # "-O0", "-g",
            ],
            include_dirs=include_dirs,
            define_macros=define_macros,
            library_dirs=library_dirs,
            libraries=libraries,
            language='c++',
            runtime_library_dirs=[],  # CUDA['lib64']],
    )

    setup(
            name=name,
            ext_modules=cythonize(
                    [ext],
                    gdb_debug=True,
                    language_level=3, ),
            # cmdclass={'build_ext': custom_build_ext},
            # zip_safe=False
    )


if __name__ == '__main__':
    try:
        # uncomment this line to disable the multithreaded compilation
        # import step_by_step

        from multiprocessing import Pool
        pool = Pool(maxtasksperchild=1)  # process per core
        # proces data_inputs iterable with pool
        pool.map(dependencies_module, listMod)
        # proces data_inputs iterable with pool
        pool.map(compile_module, listMod)
    except ImportError:
        for name in listMod:
            dependencies_module(name)
            compile_module(name)
'''
    finally:
        compile_module("shesha")

ext = Extension(
        'shesha',
        sources=['src/shesha.pyx'],
        library_dirs=library_dirs,
        libraries=libraries,
        language='c++',
        runtime_library_dirs=[],  # CUDA['lib64']],
        extra_compile_args={'g++': []},
        include_dirs=include_dirs, )
'''

# def locate_MPI():
#     """return mpi config.
#
#     NEED ENVIRONMENT VARIABLE: MPICXX
#
#     mpi_congig['mpicxx']   :compiler (mpi)
#     mpi_config['cincdirs'] :include directories (mpi)
#     mpi_config['clibs']    :libraries (mpi)
#     mpi_config['pincdirs'] :include directories (mpi4py)
#     """
#     import mpi4py
#
#     # add env variable to module
#     # mpicxx=os.environ["MPICXX"]
#     mpicxx = 'mpicxx'
#
#     # add env variable to module
#     cincdirs = subprocess.check_output([mpicxx, "--showme:incdirs"])
#     cincdirs = cincdirs.split()
#     clibs = subprocess.check_output([mpicxx, "--showme:libs"])
#     clibs = clibs.split()
#     clibdir = subprocess.check_output([mpicxx, "--showme:libdirs"])
#     clibdirs = []
#     clibdirs = clibdir.split()
#
#     pincdirs = [mpi4py.get_include(), numpy.get_include()]
#
#     mpi_config = {
#             'mpicxx': mpicxx, 'cincdirs': cincdirs, 'pincdirs': pincdirs,
#             'clibs': clibs, 'clibdirs': clibdirs}
#
#     return mpi_config
