#!/bin/bash
LOCAL_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

if [ -z $1 ]
then
    if [ ! -z $COMPASS_INSTALL_ROOT ]
    then
        echo "COMPASS_INSTALL_ROOT found."
        COMPASS_INSTALL_PATH=$COMPASS_INSTALL_ROOT
    else
        echo "Installation destination not found. Use default."
        COMPASS_INSTALL_PATH=$LOCAL_DIR/local
    fi
elif [[ $1 = "cleanup" ]]
then
    rm -rf $LOCAL_DIR/build $LOCAL_DIR/local
    COMPASS_INSTALL_PATH=$LOCAL_DIR/local
else
    echo "Installation destination provided."
    COMPASS_INSTALL_PATH=$1
fi

echo "Installing compass at ${COMPASS_INSTALL_PATH}"

if [ ! -z "$CUDA_SM" ]
then
    echo "CUDA_SM found: ${CUDA_SM}"
    CUDA_ARG='-o cuda_sm='${CUDA_SM}
else
    echo "CUDA_SM not found. Use auto detection."
fi

if [[ "${COMPASS_DO_HALF,,}" =~ ^(yes|true|on)$ ]]
then
    COMPASS_DO_HALF=True
else
    COMPASS_DO_HALF=False
fi

# If `conanlocal.py` exists, use it instead of `conanfile.py`
if [ -f $LOCAL_DIR/conanlocal.py ]
then
    CONAN_LOCATION=$LOCAL_DIR/conanlocal.py
else
    CONAN_LOCATION=$LOCAL_DIR
fi

# check if need to be built
BUILD_LIBS="${BUILD_LIBS:=True}"

if [[ -z $PYTHON_VERSION ]]
then
    PYTHON_VERSION=$(python --version | cut -d' ' -f2 | cut -d'.' -f1,2)
    echo "python version ${PYTHON_VERSION} used"
else
    echo "python version ${PYTHON_VERSION} provided"
fi

# Exit on first error.
set -e

# DEBUG="-s build_type=Debug"

# Resolves dependencies.
conan build $CONAN_LOCATION -b missing   \
    ${CUDA_ARG} ${DEBUG}                 \
    -o half=${COMPASS_DO_HALF}           \
    -o libs=${BUILD_LIBS}                \
    -o python_version=${PYTHON_VERSION}

if [[ -z $DEBUG ]]
then
    CONAN_PRESET="conan-release"
else
    CONAN_PRESET="conan-debug"
fi

cmake -DCMAKE_INSTALL_PREFIX=${COMPASS_INSTALL_PATH} --preset ${CONAN_PRESET}
cmake --build -t install --preset ${CONAN_PRESET}

echo
echo "Configuration and installation done, next time, you can simply run this command to install compass:"
echo "  cmake --build -t install --preset ${CONAN_PRESET}"
echo

# The commands generate build system in build/ subfolder, build it and install it
# conan build $CONAN_LOCATION -bf $LOCAL_DIR/build -pf ${COMPASS_INSTALL_PATH} && cmake --install build

# Only use conan to configure due to no paralelism buring build (bug).
# conan build $CONAN_LOCATION -bf $LOCAL_DIR/build -pf ${COMPASS_INSTALL_PATH} --configure && cmake --build build -j

# Install compass to specified location
# conan package $CONAN_LOCATION -bf $LOCAL_DIR/build -pf ${COMPASS_INSTALL_PATH}
