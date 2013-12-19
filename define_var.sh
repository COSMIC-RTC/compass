
# CUDA default definitions
#${STATE?"Need to set STATE"}
export CUDA_ROOT=/usr/local/cuda
export CUDA_INC_PATH=$CUDA_ROOT/include
export CUDA_LIB_PATH=$CUDA_ROOT/lib
export CUDA_LIB_PATH_64=$CUDA_ROOT/lib64
export CPLUS_INCLUDE_PATH=$CUDA_INC_PATH

# CULA default definitions
export CULA_ROOT=/usr/local/cula
export CULA_INC_PATH=$CULA_ROOT/include
export CULA_LIB_PATH=$CULA_ROOT/lib
export CULA_LIB_PATH_64=$CULA_ROOT/lib64

# YOGA default definitions
export COMPASS_ROOT_DIR=$(pwd)
export YOGA_DIR=$COMPASS_ROOT_DIR/yoga_ao
export YOGA_AO_DIR=$YOGA_DIR/yoga_ao
export YOGA_AO_TOP=$YOGA_AO_DIR

export LD_LIBRARY_PATH=$COMPASS_ROOT_DIR/libcarma:$COMPASS_ROOT_DIR/libsutra:$CUDA_LIB_PATH_64:$CUDA_LIB_PATH:$CULA_LIB_PATH_64:$CULA_LIB_PATH:$LD_LIBRARY_PATH

export YORICK_PATH=~/yorick.git/relocate/bin
export PATH=$YORICK_PATH:$PATH
#alias yorick="rlwrap yorick"

echo "COMPASS will use this configuration, please add those lines into the .bashrc"

echo "# CUDA default definitions"
echo "export CUDA_ROOT=$CUDA_ROOT"
echo "export CUDA_INC_PATH=\$CUDA_ROOT/include"
echo "export CUDA_LIB_PATH=\$CUDA_ROOT/lib"
echo "export CUDA_LIB_PATH_64=\$CUDA_ROOT/lib64"
echo "export CPLUS_INCLUDE_PATH=\$CUDA_INC_PATH"
echo
echo "# CULA default definitions"
echo "export CULA_ROOT=$CULA_ROOT"
echo "export CULA_INC_PATH=\$CULA_ROOT/include"
echo "export CULA_LIB_PATH=\$CULA_ROOT/lib"
echo "export CULA_LIB_PATH_64=\$CULA_ROOT/lib64"
echo
echo "# YOGA default definitions"
echo "export COMPASS_ROOT_DIR=$COMPASS_ROOT_DIR"
echo "export YOGA_DIR=\$COMPASS_ROOT_DIR/yoga_ao"
echo "export YOGA_AO_DIR=\$YOGA_DIR/yoga_ao"
echo "export YOGA_AO_TOP=\$YOGA_AO_DIR"
echo
echo "export LD_LIBRARY_PATH=\$COMPASS_ROOT_DIR/libcarma:\$COMPASS_ROOT_DIR/libsutra:\$CUDA_LIB_PATH_64:\$CUDA_LIB_PATH:\$CULA_LIB_PATH_64:\$CULA_LIB_PATH:\$LD_LIBRARY_PATH"
echo
echo "export YORICK_PATH=$YORICK_PATH"
echo "export PATH=\$YORICK_PATH:\$PATH"
echo 'alias yorick="rlwrap yorick"'
