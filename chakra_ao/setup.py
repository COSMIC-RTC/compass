import sys
print sys.prefix
print "SYS.PATH"
print "======================================"
print sys.path
print "======================================"

import  os, re
from os.path import join as pjoin
import subprocess
from distutils.core import setup
#from Cython.Build import cythonize
from Cython.Distutils import build_ext

## enable to dump annotate html for each pyx source file
import Cython.Compiler.Options
#Cython.Compiler.Options.annotate = True


import numpy
from distutils.extension import Extension

try:
    import mpi4py
    MPI4PY=1
except ImportError:
    MPI4PY=0


chakra_path=os.environ.get('CHAKRA')
if(chakra_path is None):
    raise EnvironmentError("Environment variable 'CHAKRA' must be define")
sys.path.append(chakra_path+'/src')



chakra_ao_path=os.environ.get('CHAKRA_AO')
if(chakra_ao_path is None):
    raise EnvironmentError("Environment variable 'CHAKRA_AO' must be define")

def find_in_path(name, path):
    "Find a file in a search path"
    #adapted fom http://code.activestate.com/recipes/52224-find-a-file-given-a-search-path/
    for dir in path.split(os.pathsep):
        binpath = pjoin(dir, name)
        if os.path.exists(binpath):
            return os.path.abspath(binpath)
    return None



def locate_cuda():
    """Locate the CUDA environment on the system

    Returns a dict with keys 'home', 'nvcc', 'include', and 'lib64'
    and values giving the absolute path to each directory.

    Starts by looking for the CUDAHOME env variable. If not found, everything
    is based on finding 'nvcc' in the PATH.
    """

    # first check if the CUDA_ROOT env variable is in use
    if 'CUDA_ROOT' in os.environ:
        home = os.environ['CUDA_ROOT']
        nvcc = pjoin(home, 'bin', 'nvcc')
    else:
        # otherwise, search the PATH for NVCC
        nvcc = find_in_path('nvcc', os.environ['PATH'])
        if nvcc is None:
            raise EnvironmentError('The nvcc binary could not be '
                'located in your $PATH. Either add it to your path, or set $CUDAHOME')
        home = os.path.dirname(os.path.dirname(nvcc))

    cudaconfig = {'home':home, 
                  'nvcc':nvcc,
                  'include': pjoin(home, 'include'),
                  'lib64': pjoin(home, 'lib64')}
    for k, v in cudaconfig.iteritems():
        if not os.path.exists(v):
            raise EnvironmentError('The CUDA %s path could not be located in %s' % (k, v))

    return cudaconfig


def locate_compass():
    """Locate compass library
    """

    if  'COMPASS_ROOT_DIR' in os.environ:
        root_compass = os.environ['COMPASS_ROOT_DIR']
        
    else:
        raise EnvironmentError('You must load compass environment')
        
    compass_config = {'inc_sutra':root_compass+'/libsutra/include.h','inc_carma':root_compass+'/libcarma/include.h',
                      'inc_chakra':root_compass+'/chakra', 'lib':root_compass}
    
    return compass_config


CUDA=locate_cuda()
COMPASS=locate_compass()
# Obtain the numpy include directory.  This logic works across numpy versions.
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()


def customize_compiler_for_nvcc(self):
    """inject deep into distutils to customize how the dispatch
    to gcc/nvcc works.
    
    If you subclass UnixCCompiler, it's not trivial to get your subclass
    injected in, and still have the right customizations (i.e.
    distutils.sysconfig.customize_compiler) run on it. So instead of going
    the OO route, I have this. Note, it's kindof like a wierd functional
    subclassing going on."""
    
    # tell the compiler it can processes .cu
    self.src_extensions.append('.cu')

    # save references to the default compiler_so and _comple methods
    default_compiler_so = self.compiler_so
    super = self._compile

    # now redefine the _compile method. This gets executed for each
    # object but distutils doesn't have the ability to change compilers
    # based on source extension: we add it.
    def _compile(obj, src, ext, cc_args, extra_postargs, pp_opts):
        if os.path.splitext(src)[1] == '.cu':
            # use the cuda for .cu files
            self.set_executable('compiler_so', CUDA['nvcc'])
            # use only a subset of the extra_postargs, which are 1-1 translated
            # from the extra_compile_args in the Extension class
            postargs = extra_postargs['nvcc']
        else:
            postargs = extra_postargs['g++']

        super(obj, src, ext, cc_args, postargs, pp_opts)
        # reset the default compiler_so, which we might have changed for cuda
        self.compiler_so = default_compiler_so

    # inject our redefined _compile method into the class
    self._compile = _compile


def locate_MPI():
    """return mpi config

    NEED ENVIRONMENT VARIABLE: MPICXX

    mpi_congig['mpicxx']   :compiler (mpi)
    mpi_config['cincdirs'] :include directories (mpi)
    mpi_config['clibs']    :librairies (mpi)
    mpi_config['pincdirs'] :include directories (mpi4py)
    """

    #add env variable to module
    #mpicxx=os.environ["MPICXX"]
    mpicxx='mpicxx'

    #add env variable to module
    cincdirs=subprocess.check_output([mpicxx,"--showme:incdirs"])
    cincdirs=cincdirs.split()
    clibs=subprocess.check_output([mpicxx,"--showme:libs"])
    clibs=clibs.split()
    clibdir=subprocess.check_output([mpicxx,"--showme:libdirs"])
    clibdirs=[]
    clibdirs=clibdir.split()

    pincdirs=[mpi4py.get_include(),numpy.get_include()]

    mpi_config = {'mpicxx':mpicxx,'cincdirs':cincdirs, 'pincdirs':pincdirs, 'clibs':clibs,'clibdirs':clibdirs}

    return mpi_config


# run the customize_compiler
class custom_build_ext(build_ext):
    def build_extensions(self):
        if (os.path.exists('./chakra_ao.cpp')): os.remove('./chakra_ao.cpp')
        customize_compiler_for_nvcc(self.compiler)
        build_ext.build_extensions(self)


source=['chakra_ao']
librairies=['sutra']
include_dirs = [numpy_include, 
                CUDA['include'],
                COMPASS['inc_carma'],
                COMPASS['inc_sutra'],
                COMPASS['inc_chakra']]

library_dirs=[COMPASS['lib']+"/libsutra"]

if MPI4PY==1:
    MPI=locate_MPI()
    librairies.extend(MPI['clibs'])
    include_dirs.extend(MPI['cincdirs'])
    include_dirs.extend(MPI['pincdirs'])
    library_dirs.extend(MPI['clibdirs'])

ext=Extension('chakra_ao',
              sources=['src/chakra_ao.pyx'],
              library_dirs=library_dirs,
              libraries=librairies,
              language='c++',
              runtime_library_dirs=[],#CUDA['lib64']],
              extra_compile_args={'g++': []},
              include_dirs = include_dirs,
              )


setup(
    name="chakra_ao",
    ext_modules=[ext],
    cmdclass={'build_ext': custom_build_ext},
    zip_safe=False
)
