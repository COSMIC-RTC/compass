
from distutils.core import setup, Extension
import sys 

def search_on_path(filenames):
    """Find file on system path."""
    # http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/52224

    from os.path import exists, join, abspath
    from os import pathsep, environ

    search_path = environ["PATH"]

    paths = search_path.split(pathsep)
    for path in paths:
        for filename in filenames:
            if exists(join(path, filename)):
                return abspath(join(path, filename))

#nvcc_path = search_on_path(["nvcc", "nvcc.exe"])
#if nvcc_path is None:
#    print("*** CUDA_ROOT not set, and nvcc not in path. Giving up.")
#    sys.exit(1)
#

from os import pathsep, environ
module_yoga = Extension('Yoga',
                    define_macros = [('MAJOR_VERSION', '1'),
                                     ('MINOR_VERSION', '0')],
                    include_dirs = ['/usr/local/include', 
                                    '/usr/local/cuda/include', 
                                    '/usr/local/cula/include', 
                                    '%s/include.h'%(environ["LIBYOGAPATH"])],
                    extra_link_args = ['-Wl,-no-undefined',
                                       '-Wl,--warn-once',
				       '-fPIC'],
                    libraries = ['yoga', 'python%d.%d'%(sys.version_info[0],sys.version_info[1]),
                                 'stdc++', 'm'],
                    library_dirs = [environ["LIBYOGAPATH"]],
                    sources = ['yoga_py.cpp'])

setup (name = 'Yoga',
       version = '1.0',
       description = 'This is a Yoga package',
       author = 'Arnaud Sevin & Damien Gratadour',
       author_email = '',
       url = 'http://docs.python.org/extending/building',
       long_description = '''
This is really just a Yoga package.
''',
       ext_modules = [module_yoga])

