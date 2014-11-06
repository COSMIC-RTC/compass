import os
import shutil
import glob

"""
kind of make clean
 because 'python setup.py clean' does not clean everything
 remove build directory and other generated files. 
"""

for dirname in ('build', '__pycache__'):
    if os.path.exists(dirname):
        shutil.rmtree(dirname)
        
for filepattern in ('*.so', '*.pyc'):
    for filename in glob.glob(filepattern):
        
            print("delete {}".format(filename))        
            os.remove(filename)

if os.path.exists('wrapper_chakra.cpp'):
    os.remove('wrapper_chakra.cpp')
if os.path.exists('wrapper_chakra.pyx'):
    os.remove('wrapper_chakra.pyx')
