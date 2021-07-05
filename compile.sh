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

if [ ! -z $CUDA_SM ]
then
    echo "CUDA_SM found: ${CUDA_SM}"
else
    echo "CUDA_SM not found. Use auto detection."
    CUDA_SM=Auto
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

if [ -z $PYTHON_VERSION ]
then
    PYTHON_VERSION=$(python --version | cut -d' ' -f2 | cut -d'.' -f1,2)
    echo "python version ${PYTHON_VERSION} used"
fi

# Exit on first error.
set -e

# Resolves dependencies.
conan install $CONAN_LOCATION -if $LOCAL_DIR/build -b missing \
    -o emu:cuda_sm=${CUDA_SM}                                 \
    -o compass:half=${COMPASS_DO_HALF}                        \
    -o compass:python_version=${PYTHON_VERSION}

# Build compass
# Only use conan to configure due to no paralelism buring build (bug).
conan build $CONAN_LOCATION -bf $LOCAL_DIR/build -pf ${COMPASS_INSTALL_PATH} --configure
cmake --build build -j
# Install compass to specified location
conan package $CONAN_LOCATION -bf $LOCAL_DIR/build -pf ${COMPASS_INSTALL_PATH}
