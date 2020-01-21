#!/bin/bash

if [ -z $1 ]
then
    PYTHON_INSTALL_PATH=$COMPASS_ROOT/local
else
    PYTHON_INSTALL_PATH=$1
fi

if [ ! -d build ]
then
        echo "Create build directory"
        mkdir -p build
fi

if [ -z $COMPASS_DO_HALF ]
then
    COMPASS_DO_HALF="OFF"
fi

cd build

# BUILD_TOOL="-GNinja" # build with ninja instead of make
# COMPASS_DEBUG="-DCMAKE_BUILD_TYPE=Debug"
NCPUS=`fgrep processor /proc/cpuinfo | wc -l`

cmake .. -DCMAKE_INSTALL_PREFIX=$PYTHON_INSTALL_PATH -Ddo_half=$COMPASS_DO_HALF $COMPASS_DEBUG $BUILD_TOOL
cmake --build . --target install -- -j $NCPUS
