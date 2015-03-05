import  os, re
from os.path import join as pjoin
from distutils.core import setup
#from Cython.Build import cythonize
from Cython.Distutils import build_ext

import numpy
from distutils.extension import Extension



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
        inc_compass  = pjoin(root_compass,)#os.environ['COMPASS_ROOT_DIR']
        lib_compass  = pjoin(root_compass,)#os.environ['COMPASS_ROOT_DIR']
        
    else:
        raise EnvironmentError('You must load compass environment')
        
    compass_config = {'inc':inc_compass, 'lib':lib_compass}
    
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
            postargs = extra_postargs['gcc']

        super(obj, src, ext, cc_args, postargs, pp_opts)
        # reset the default compiler_so, which we might have changed for cuda
        self.compiler_so = default_compiler_so

    # inject our redefined _compile method into the class
    self._compile = _compile


# run the customize_compiler
class custom_build_ext(build_ext):
    def build_extensions(self):
        if (os.path.exists('./chakra_ao.cpp')): os.remove('./chakra_ao.cpp')
        customize_compiler_for_nvcc(self.compiler)
        build_ext.build_extensions(self)



ext=Extension('chakra_ao',
              sources=['chakra_ao.pyx'],
              library_dirs=[CUDA['lib64'],COMPASS['lib']+"/libsutra"],#,COMPASS['lib']+"/libcarma"],
              libraries=[ 'sutra'],
              language='c++',
              runtime_library_dirs=[],#CUDA['lib64']],
              extra_compile_args={'gcc': []},
              include_dirs = [numpy_include, 
                              CUDA['include'],
                              COMPASS['inc']+'/libcarma/include.h',
                              COMPASS['inc']+'/libsutra/include.h',
                              '../chakra'])

setup(
    name="chakra_ao",
    ext_modules=[ext],
    cmdclass={'build_ext': custom_build_ext},
    zip_safe=False
)
