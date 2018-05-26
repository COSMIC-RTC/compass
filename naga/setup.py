import os
import re
from os.path import join as pjoin
from distutils.core import setup, Extension
import numpy
import sys
from distutils.command.clean import clean as _clean
from distutils.dir_util import remove_tree

from Cython.Build import cythonize
# import shutil

# Remove the "-Wstrict-prototypes" compiler option, which isn't valid for C++.
import distutils.sysconfig
cfg_vars = distutils.sysconfig.get_config_vars()
for key, value in list(cfg_vars.items()):
    if type(value) == str:
        cfg_vars[key] = value.replace("-Wstrict-prototypes", "")
# ==================================

listMod = ['context', 'streams', 'host_obj', 'obj', 'magma', 'timer', 'sparse_obj']
dependencies = {
        'streams': ['context'],
        'host_obj': ['context', 'streams'],
        'obj': ['context', 'streams'],
        'magma': ['obj', 'host_obj'],
        'timer': ['context'],
        'sparse_obj': ['context', 'obj']
}

naga_path = os.environ.get('NAGA_ROOT')
if (naga_path is None):
    raise EnvironmentError("Environment variable 'NAGA_ROOT' must be define")

# deprecated
# def find_in_path(name, path):
#     "Find a file in a search path"
#     # adapted fom http://code.activestate.com/recipes/52224-find-a-file-given-a-search-path/
#     for dir in path.split(os.pathsep):
#         binpath = pjoin(dir, name)
#         if os.path.exists(binpath):
#             return os.path.abspath(binpath)
#     return None


def locate_compass():
    """Locate compass library
    """

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

# --------------------------------------------------------------------
# Clean target redefinition - force clean everything
# --------------------------------------------------------------------
relist = ['^.*~$', '^core\.*$', '^#.*#$', '^.*\.aux$', '^.*\.pyc$', '^.*\.o$']
reclean = []

for restring in relist:
    reclean.append(re.compile(restring))


def wselect(args, dirname, names):
    for n in names:
        for rev in reclean:
            if (rev.match(n)):
                os.remove("%s/%s" % (dirname, n))
                break


class clean(_clean):

    def walkAndClean(self):
        os.path.walk("..", wselect, [])

    def run(self):
        module_lib = pjoin('.', 'naga.so')
        if (os.path.exists(module_lib)):
            os.remove(module_lib)
        if (os.path.exists('./naga')):
            os.remove('./naga')
        if (os.path.exists('./naga/*.cpp')):
            os.remove('./naga/*.cpp')
        if (os.path.exists('./naga/obj.pyx')):
            os.remove('./naga/obj.pyx')
        if (os.path.exists('./naga/host_obj.pyx')):
            os.remove('./naga/host_obj.pyx')
        if (os.path.exists('./naga/magma.pyx')):
            os.remove('./naga/magma.pyx')
        if (os.path.exists('./build')):
            remove_tree('./build')
        if (os.path.exists('./dist')):
            remove_tree('./dist')
        self.walkAndClean()


library_dirs = [COMPASS['lib'] + '/libcarma']
libraries = ['carma']

include_dirs = [numpy_include, COMPASS['inc_carma']]

if 'MKLROOT' in os.environ:
    mkl_root = os.environ.get('MKLROOT')
    print("mkl_root: ", mkl_root)
    library_dirs.append(mkl_root + '/mkl/lib/intel64/')
    libraries.append('mkl_mc3')
    libraries.append('mkl_def')

if 'CUDA_INC_PATH' in os.environ:
    include_dirs.append(os.environ['CUDA_INC_PATH'])
    library_dirs.append(os.environ['CUDA_LIB_PATH_64'])
    libraries.append('cudart')
else:
    raise EnvironmentError("Environment variable 'CUDA_INC_PATH' must be define")

# print("library_dirs", library_dirs)
# print("libraries", libraries)

#######################
#  extension
#######################
# deprecated
# ext = Extension('naga',
#                sources=['src/naga.pyx'
#                    ],
#                library_dirs=library_dirs,
#                libraries=libraries,
#                language='c++',
#                runtime_library_dirs=[],
#                # this syntax is specific to this build system
#                # we're only going to use certain compiler args with nvcc and not with gcc
#                # the implementation of this trick is in customize_compiler() below
#                #extra_compile_args=["-O0", "-g"],
#                #extra_compile_args={'g++': [],},
#                                    #nvcc not needed (cuda code alreay compiled)
#                                    #'nvcc': ['-gencode '+os.environ['GENCODE'],
#                                    #         '--ptxas-options=-v',
#                                    #         '-c', '--compiler-options',
#                                    #         "'-fPIC'"]},
#                include_dirs = [numpy_include,
#                                COMPASS['inc']+'/libcarma/include.h'])
#

# deprecated
# def customize_compiler_for_nvcc(self):
#     """inject deep into distutils to customize how the dispatch
#     to gcc/nvcc works.
#
#     If you subclass UnixCCompiler, it's not trivial to get your subclass
#     injected in, and still have the right customizations (i.e.
#     distutils.sysconfig.customize_compiler) run on it. So instead of going
#     the OO route, I have this. Note, it's kindof like a wierd functional
#     subclassing going on."""
#
#     # save references to the default compiler_so and _comple methods
#     default_compiler_so = self.compiler_so
#     super = self._compile
#
#     # now redefine the _compile method. This gets executed for each
#     # object but distutils doesn't have the ability to change compilers
#     # based on source extension: we add it.
#     def _compile(obj, src, ext, cc_args, extra_postargs, pp_opts):
#         postargs = extra_postargs['g++']
#
#         super(obj, src, ext, cc_args, postargs, pp_opts)
#         # reset the default compiler_so, which we might have changed for cuda
#         self.compiler_so = default_compiler_so
#
#     # inject our redefined _compile method into the class
#     self._compile = _compile

# deprecated
# run the customize_compiler
# class custom_build_ext(build_ext):
#     def build_extensions(self):
#         if (os.path.exists('./naga.cpp')):
#             os.remove('./naga.cpp')
#         customize_compiler_for_nvcc(self.compiler)
#         build_ext.build_extensions(self)

# dal with generated sources files
if 'build_ext' in sys.argv or 'develop' in sys.argv or 'install' in sys.argv:
    generator = os.path.join(os.path.abspath('.'), 'naga/process_tmpl.py')
    d = {'__file__': generator}
    exec(compile(open(generator).read(), generator, 'exec'), d)
    d['main'](None)


def dependencies_module(name):
    print("=======================================")
    print("resolving dependencies for", name)
    print("=======================================")
    try:
        dep = dependencies[name]
        print("dependencies:", dep)
        if (os.path.exists("naga/" + name + ".cpp")):
            for d in dep:
                if (os.stat("naga/" + d + ".pyx").st_mtime >
                            os.stat("naga/" + name + ".cpp").st_mtime):
                    # cpp file outdated
                    # cpp file outdated if exists
                    if (os.path.exists("naga/" + name + ".cpp")):
                        os.remove("naga/" + name + ".cpp")
    except KeyError as e:
        print(e)


def compile_module(name):
    print("=======================================")
    print("creating module ", name)
    print("=======================================")
    ext = Extension(
            "naga." + name,
            sources=[naga_path + '/naga/' + name + '.pyx'],
            extra_compile_args=[
                    "-Wno-unused-function",
                    "-Wno-unused-label",
                    "-Wno-cpp",
                    "-std=c++11",
                    # "-O0", "-g",
            ],
            include_dirs=include_dirs,
            library_dirs=library_dirs,
            libraries=libraries,
            language='c++',
            runtime_library_dirs=[],  # CUDA['lib64']],
            define_macros=[], )

    setup(
            name=name,
            ext_modules=cythonize(
                    [ext],
                    # gdb_debug=True,
                    language_level=3, )
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
    finally:
        compile_module("naga")

# setup(name='naga',
#
#      ext_modules = cythonize([ext]),
#
#      # inject our custom trigger
#      #cmdclass={'build_ext': custom_build_ext},
#
#      # since the package has c code, the egg cannot be zipped
#      zip_safe=False)
