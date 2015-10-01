import  os, re
from os.path import join as pjoin
from distutils.core import setup
#from setuptools import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy
import sys
from distutils.command.clean import clean as _clean
from distutils.dir_util import remove_tree

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



CUDA = locate_cuda()
COMPASS = locate_compass()


# Obtain the numpy include directory.  This logic works across numpy versions.
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()


# --------------------------------------------------------------------
# Clean target redefinition - force clean everything
# --------------------------------------------------------------------
relist=['^.*~$','^core\.*$','^#.*#$','^.*\.aux$','^.*\.pyc$','^.*\.o$']
reclean=[]

for restring in relist:
  reclean.append(re.compile(restring))

def wselect(args,dirname,names):
  for n in names:
    for rev in reclean:
      if (rev.match(n)):
        os.remove("%s/%s"%(dirname,n))
        break

class clean(_clean):
  def walkAndClean(self):
    os.path.walk("..",wselect,[])
  def run(self):
    module_lib = pjoin('.',module_name+'.so')
    if (os.path.exists(module_lib)): os.remove(module_lib)
    if (os.path.exists('./src/wrapper_chakra.cpp')): os.remove('./src/wrapper_chakra.cpp')
    if (os.path.exists('./src/wrapper_chakra_obj.pyx')): os.remove('./src/wrapper_chakra_obj.pyx')
    if (os.path.exists('./src/wrapper_chakra_host_obj.pyx')): os.remove('./src/wrapper_chakra_host_obj.pyx')
    if (os.path.exists('./src/wrapper_magma.pyx')): os.remove('./src/wrapper_magma.pyx')
    if (os.path.exists('./build')): remove_tree('./build')
    if (os.path.exists('./dist')):  remove_tree('./dist')
    self.walkAndClean()

#######################
#  extension
#######################
ext = Extension('chakra',
                sources=[#COMPASS['inc']+'/libcarma/src.cpp/carma_context.cpp',
                    #COMPASS['inc']+'/libcarma/src.cpp/carma_cublas.cpp',
                    #COMPASS['inc']+'/libcarma/src.cpp/carma_fft_conv.cpp',
                    #'wrapper_chakra_obj.pyx',
                    #'wrapper_chakra_host_obj.pyx'
                    'src/wrapper_chakra.pyx'
                    ],
                library_dirs=[CUDA['lib64'], COMPASS['lib']+"/libcarma", "/usr/local/intel/composer_xe_2013_sp1.1.106/mkl"],
                libraries=['cudart', 'cufft','cublas','carma',
					"mkl_mc3","mkl_def"],
                language='c++',
                runtime_library_dirs=[CUDA['lib64']],
                # this syntax is specific to this build system
                # we're only going to use certain compiler args with nvcc and not with gcc
                # the implementation of this trick is in customize_compiler() below
                extra_compile_args={'gcc': [],},
                                    #nvcc not needed (cuda code alreay compiled)
                                    #'nvcc': ['-gencode '+os.environ['GENCODE'], 
                                    #         '--ptxas-options=-v', 
                                    #         '-c', '--compiler-options', 
                                    #         "'-fPIC'"]},
                include_dirs = [numpy_include, 
                                CUDA['include'], 
                                COMPASS['inc']+'/libcarma/include.h'])



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
        if (os.path.exists('./wrapper_chakra.cpp')): os.remove('./wrapper_chakra.cpp')
        customize_compiler_for_nvcc(self.compiler)
        build_ext.build_extensions(self)


# dal with generated sources files
if 'build_ext' in sys.argv or 'develop' in sys.argv or 'install' in sys.argv:
    generator = os.path.join( os.path.abspath('.'), 'src/process_tmpl.py')
    d = {'__file__': generator }
    execfile(generator, d)
    d['main'](None)
                              



setup(name='chakra',

      ext_modules = [ext],

      # inject our custom trigger
      cmdclass={'build_ext': custom_build_ext},

      # since the package has c code, the egg cannot be zipped
      zip_safe=False)
