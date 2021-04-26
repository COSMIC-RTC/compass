#!/bin/bash
LOCAL_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

if [ ! -z $1 ]
then
    COMPASS_INSTALL_ROOT=$1
fi

if [ -z $COMPASS_INSTALL_ROOT ]
then
    COMPASS_INSTALL_ROOT=$LOCAL_DIR/local
fi

if [ -z $COMPASS_DO_HALF ]
then
    COMPASS_DO_HALF="OFF"
fi

# If `conanlocal.txt` exists, use it instead of `conanfile.txt`
if [ -f $LOCAL_DIR/conanlocal.py ]
then
    CONAN_LOCATION=conanlocal.py
else
    CONAN_LOCATION=$LOCAL_DIR
fi

# BUILD_TOOL="-GNinja" # build with ninja instead of make
# COMPASS_DEBUG="-DCMAKE_BUILD_TYPE=Debug"
#NCPUS=`fgrep processor /proc/cpuinfo | wc -l`

conan install -if build --build=missing $CONAN_LOCATION
# conan build -bf build  $CONAN_LOCATION
# conan package -bf build -pf $COMPASS_INSTALL_ROOT $CONAN_LOCATION

cd build
cmake $LOCAL_DIR -DPYTHON_EXECUTABLE=$(which python) -DCMAKE_INSTALL_PREFIX=$COMPASS_INSTALL_ROOT -Ddo_half=$COMPASS_DO_HALF $COMPASS_DEBUG $BUILD_TOOL
cmake --build . --target install -- -j $NCPUS

# conan export-pkg . conan/stable -f
# conan upload compass/5.0@conan/stable --all -r=hippo6
