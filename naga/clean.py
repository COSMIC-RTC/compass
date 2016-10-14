import os
import shutil
import glob

"""
kind of make clean
 because 'python setup.py clean' does not clean everything
 remove build directory and other generated files. 
"""

SRC='src/'

for dirname in ('build', '__pycache__'):
    if os.path.exists(dirname):
        shutil.rmtree(dirname)
        
for filepattern in ('*.so', SRC+'*.pyc'):
    for filename in glob.glob(filepattern):
        
            print("delete {}".format(filename))        
            os.remove(filename)

FILE=['naga_obj', 'naga_host_obj', 'naga_sparse_obj', 'naga_magma']

if os.path.exists(SRC+'naga.cpp'):
    os.remove(SRC+'naga.cpp')
for f in FILE:
    if (os.path.exists(SRC+f+'.pyx') and os.path.exists(SRC+f+'.pyx.in') ) :
        os.remove(SRC+f+'.pyx')
    if os.path.exists(SRC+f+'.cpp'):
        os.remove(SRC+f+'.cpp')
