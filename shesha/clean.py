"""Kind of make clean.

because 'python setup.py clean' does not clean everything
remove build directory and other generated files.
"""

import os
import shutil
import glob

SRC = 'shesha/sutra_bind'

for dirname in ('build', '__pycache__'):
    if os.path.exists(dirname):
        shutil.rmtree(dirname)

for filepattern in ('*.so', SRC + '/*.so', SRC + '/*.pyc', SRC + '/*.cpp'):
    for filename in glob.glob(filepattern):

        print(("delete {}".format(filename)))
        os.remove(filename)

if os.path.exists(SRC + 'shesha.cpp'):
    os.remove(SRC + 'shesha.cpp')
