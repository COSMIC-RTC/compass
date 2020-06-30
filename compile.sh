#!/bin/bash
LOCAL_DIR="$(realpath --relative-to=$(pwd) $( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd ))"

if [ -z $1 ]
then
    PYTHON_INSTALL_PATH=$COMPASS_ROOT/local
else
    PYTHON_INSTALL_PATH=$1
fi

if [ ! -d $LOCAL_DIR/build ]
then
    echo "Create build directory"
    mkdir -p $LOCAL_DIR/build
fi

if [ -z $COMPASS_DO_HALF ]
then
    COMPASS_DO_HALF="OFF"
fi

# If `conanlocal.txt` exists, use it instead of `conanfile.txt`
if [ -f $LOCAL_DIR/conanlocal.py ]
then
    CONAN_LOCATION=../conanlocal.py
else
    CONAN_LOCATION=..
fi

cd $LOCAL_DIR/build

# BUILD_TOOL="-GNinja" # build with ninja instead of make
# COMPASS_DEBUG="-DCMAKE_BUILD_TYPE=Debug"
NCPUS=`fgrep processor /proc/cpuinfo | wc -l`

conan install $CONAN_LOCATION --build=missing
cmake .. -DCMAKE_INSTALL_PREFIX=$PYTHON_INSTALL_PATH -Ddo_half=$COMPASS_DO_HALF $COMPASS_DEBUG $BUILD_TOOL
cmake --build . --target install -- -j $NCPUS
