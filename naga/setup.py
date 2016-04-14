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


listMod=['naga_context','naga_streams','naga_obj','naga_host_obj','naga_magma','naga_timer','naga_sparse_obj']
dependencies={'naga_streams':['naga_context'],
              'naga_obj':['naga_context','naga_streams'],
              'naga_host_obj':['naga_context','naga_streams'],
              'naga_magma':['naga_obj','naga_host_obj'],
              'naga_timer':[],
              'naga_sparse_obj':['naga_context','naga_obj']}

naga_path=os.environ.get('NAGA_ROOT')
if(naga_path is None):
    raise EnvironmentError("Environment variable 'NAGA_ROOT' must be define")
sys.path.append(naga_path+'/src')


def find_in_path(name, path):
    "Find a file in a search path"
    #adapted fom http://code.activestate.com/recipes/52224-find-a-file-given-a-search-path/
    for dir in path.split(os.pathsep):
        binpath = pjoin(dir, name)
        if os.path.exists(binpath):
            return os.path.abspath(binpath)
    return None

    
def locate_compass():
    """Locate compass library
    """

    if  'COMPASS_ROOT' in os.environ:
        root_compass = os.environ['COMPASS_ROOT']
        
    else:
        raise EnvironmentError("Environment variable 'COMPASS_ROOT' must be define")
        
    compass_config = {'inc_sutra':root_compass+'/libsutra/include.h','inc_carma':root_compass+'/libcarma/include.h',
                      'inc_naga':root_compass+'/naga', 'lib':root_compass}
    
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
    if (os.path.exists('./src/naga.cpp')): os.remove('./src/naga.cpp')
    if (os.path.exists('./src/naga_obj.pyx')): os.remove('./src/naga_obj.pyx')
    if (os.path.exists('./src/naga_host_obj.pyx')): os.remove('./src/naga_host_obj.pyx')
    if (os.path.exists('./src/naga_magma.pyx')): os.remove('./src/naga_magma.pyx')
    if (os.path.exists('./build')): remove_tree('./build')
    if (os.path.exists('./dist')):  remove_tree('./dist')
    self.walkAndClean()


library_dirs=[COMPASS['lib']+'/libcarma']
libraries=['carma']

print "library_dirs",library_dirs
print "libraries",libraries

mkl_root=os.environ.get('MKLROOT')
print "mkl_root:"
if(mkl_root==""):
    mkl_root=None
if(mkl_root is not None):
    print mkl_root
    library_dirs.append(mkl_root+'/mkl/lib/intel64/')
    libraries.append('mkl_mc3')
    libraries.append('mkl_def')

print "library_dirs",library_dirs
print "libraries",libraries

if  'CUDA_INC_PATH' in os.environ:
    cuda_include = os.environ['CUDA_INC_PATH']
else:
    raise EnvironmentError("Environment variable 'CUDA_INC_PATH' must be define")

include_dirs=[numpy_include, 
                COMPASS['inc_carma'],
                cuda_include
            ]

#######################
#  extension
#######################
#ext = Extension('naga',
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


def customize_compiler_for_nvcc(self):
    """inject deep into distutils to customize how the dispatch
    to gcc/nvcc works.
    
    If you subclass UnixCCompiler, it's not trivial to get your subclass
    injected in, and still have the right customizations (i.e.
    distutils.sysconfig.customize_compiler) run on it. So instead of going
    the OO route, I have this. Note, it's kindof like a wierd functional
    subclassing going on."""
    
    # save references to the default compiler_so and _comple methods
    default_compiler_so = self.compiler_so
    super = self._compile

    # now redefine the _compile method. This gets executed for each
    # object but distutils doesn't have the ability to change compilers
    # based on source extension: we add it.
    def _compile(obj, src, ext, cc_args, extra_postargs, pp_opts):
        postargs = extra_postargs['g++']

        super(obj, src, ext, cc_args, postargs, pp_opts)
        # reset the default compiler_so, which we might have changed for cuda
        self.compiler_so = default_compiler_so

    # inject our redefined _compile method into the class
    self._compile = _compile


# run the customize_compiler
class custom_build_ext(build_ext):
    def build_extensions(self):
        if (os.path.exists('./naga.cpp')): os.remove('./naga.cpp')
        customize_compiler_for_nvcc(self.compiler)
        build_ext.build_extensions(self)


# dal with generated sources files
if 'build_ext' in sys.argv or 'develop' in sys.argv or 'install' in sys.argv:
    generator = os.path.join( os.path.abspath('.'), 'src/process_tmpl.py')
    d = {'__file__': generator }
    execfile(generator, d)
    d['main'](None)
                              


from Cython.Build import cythonize

def compile_module(name):
    if(os.path.exists(naga_path+"/lib/"+name+".so") and name != "naga"):
        shutil.move(naga_path+"/lib/"+name+".so",naga_path+"/"+name+".so")
    print "======================================="
    print "creating module ",name
    print "======================================="
    try:
        dep=dependencies[name]
        print "dependencies:",dep
        if(os.path.exists("src/"+name+".cpp")):
            for d in dep:
                if (os.stat("src/"+d+".pyx").st_mtime > 
                    os.stat("src/"+name+".cpp").st_mtime):
                    #cpp file outdated
                    os.remove("src/"+name+".cpp")
    except KeyError, e:
        print e


    ext=Extension(name,
                  sources=['src/'+name+'.pyx'],
                  library_dirs=library_dirs,
                  libraries=libraries,
                  language='c++',
                  runtime_library_dirs=[],#CUDA['lib64']],
                  #extra_compile_args=["-O0", "-g"],
                  #extra_compile_args={'g++': []},
                  include_dirs = include_dirs,
                  define_macros = [],
                  )


    setup(
        name=name,
        ext_modules=cythonize([ext]),
        #cmdclass={'build_ext': custom_build_ext},
        zip_safe=False
    )
    if(os.path.exists(naga_path+"/"+name+".so") and name != "naga"):
        shutil.move(naga_path+"/"+name+".so",naga_path+"/lib/"+name+".so")

if __name__ == '__main__':
    try :
        #uncomment this line to disable the multithreaded compilation
        #import step_by_step 
        
        from multiprocessing import Pool
        pool = Pool(maxtasksperchild=1) # process per core
        pool.map(compile_module, listMod)  # proces data_inputs iterable with poo
    except ImportError:
        for name in listMod:
            compile_module(name)
    finally:
        compile_module("naga")


#setup(name='naga',
#
#      ext_modules = cythonize([ext]),
#
#      # inject our custom trigger
#      #cmdclass={'build_ext': custom_build_ext},
#
#      # since the package has c code, the egg cannot be zipped
#      zip_safe=False)
