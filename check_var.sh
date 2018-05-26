#!/bin/bash

function check_var() {
    if [[ $2 != "" ]]
    then
	echo $1 defined to $2
    else
	echo $1 not defined !!!
    fi
}

echo *** check CUDA definitions
check_var CUDA_ROOT $CUDA_ROOT
check_var CUDA_INC_PATH $CUDA_INC_PATH
check_var CUDA_LIB_PATH $CUDA_LIB_PATH
check_var CUDA_LIB_PATH_64 $CUDA_LIB_PATH_64
check_var CUDA_SM $CUDA_SM

echo *** check MAGMA definitions
check_var MAGMA_ROOT $MAGMA_ROOT

# echo *** CULA default definitions
# check_var CULA_ROOT $CULA_ROOT
# check_var CULA_INC_PATH $CULA_INC_PATH
# check_var CULA_LIB_PATH $CULA_LIB_PATH
# check_var CULA_LIB_PATH_64 $CULA_LIB_PATH_64

echo *** check COMPASS definitions
check_var COMPASS_ROOT $COMPASS_ROOT
check_var NAGA_ROOT $NAGA_ROOT
check_var SHESHA_ROOT $SHESHA_ROOT

echo *** check environment definitions
check_var PATH $PATH
check_var LD_LIBRARY_PATH $LD_LIBRARY_PATH
check_var PKG_CONFIG_PATH $PKG_CONFIG_PATH
check_var PYTHONPATH $PYTHONPATH
