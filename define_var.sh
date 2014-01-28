
# CUDA default definitions
#export GENCODE='"arch=compute_12,code=sm_12" -DGPUSHMEM=120' #for old tesla
export GENCODE='"arch=compute_20,code=sm_20" -DGPUSHMEM=200' #for Fermi
#export GENCODE='"arch=compute_30,code=sm_30" -DGPUSHMEM=300' #for Kepler 
#export GENCODE='"arch=compute_35,code=sm_35" -DGPUSHMEM=300' #for K40
#${GENCODE?"Need to set GENCODE"}

export CUDA_ROOT=/usr/local/cuda
export CUDA_INC_PATH=$CUDA_ROOT/include
export CUDA_LIB_PATH=$CUDA_ROOT/lib
export CUDA_LIB_PATH_64=$CUDA_ROOT/lib64
export CPLUS_INCLUDE_PATH=$CUDA_INC_PATH
export PATH=$CUDA_ROOT/bin:$PATH

# CULA default definitions
#export CULA_ROOT=/usr/local/cula
#export CULA_INC_PATH=$CULA_ROOT/include
#export CULA_LIB_PATH=$CULA_ROOT/lib
#export CULA_LIB_PATH_64=$CULA_ROOT/lib64

# MAGMA definitions (uncomment this line if MAGMA is installed)
export MAGMA_ROOT=/usr/local/magma
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MAGMA_ROOT/lib

# YOGA default definitions
export COMPASS_ROOT_DIR=$(pwd)
export YOGA_DIR=$COMPASS_ROOT_DIR/yoga
export YOGA_AO_DIR=$COMPASS_ROOT_DIR/yoga_ao
export YOGA_AO_TOP=$YOGA_AO_DIR

export LD_LIBRARY_PATH=$COMPASS_ROOT_DIR/libcarma:$COMPASS_ROOT_DIR/libsutra:$CUDA_LIB_PATH_64:$CUDA_LIB_PATH:$CULA_LIB_PATH_64:$CULA_LIB_PATH:$LD_LIBRARY_PATH

export YORICK_PATH=~/yorick.git/relocate/bin
export PATH=$YORICK_PATH:$PATH
#alias yorick="rlwrap -s 2000 yorick"

echo "COMPASS will use this configuration, please add those lines into the .bashrc"

echo "# CUDA default definitions"
echo "export GENCODE=\"$GENCODE\""
echo "# The value of GENCODE should be modified according to your GPU generation"
echo "export CUDA_ROOT=$CUDA_ROOT"
echo "export CUDA_INC_PATH=\$CUDA_ROOT/include"
echo "export CUDA_LIB_PATH=\$CUDA_ROOT/lib"
echo "export CUDA_LIB_PATH_64=\$CUDA_ROOT/lib64"
echo "export CPLUS_INCLUDE_PATH=\$CUDA_INC_PATH"
echo "export PATH=\$CUDA_ROOT/bin:\$PATH"
echo
echo "# CULA default definitions (optionnal)"
echo "export CULA_ROOT=$CULA_ROOT"
echo "export CULA_INC_PATH=\$CULA_ROOT/include"
echo "export CULA_LIB_PATH=\$CULA_ROOT/lib"
echo "export CULA_LIB_PATH_64=\$CULA_ROOT/lib64"
echo
echo "# MAGMA definitions (optionnal)"
echo "export MAGMA_ROOT=$MAGMA_ROOT"
echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:\$MAGMA_ROOT/lib"
echo 
echo "# YOGA default definitions"
echo "export COMPASS_ROOT_DIR=$COMPASS_ROOT_DIR"
echo "export YOGA_DIR=\$COMPASS_ROOT_DIR/yoga"
echo "export YOGA_AO_DIR=\$COMPASS_ROOT_DIR/yoga_ao"
echo "export YOGA_AO_TOP=\$YOGA_AO_DIR"
echo
echo "export LD_LIBRARY_PATH=\$COMPASS_ROOT_DIR/libcarma:\$COMPASS_ROOT_DIR/libsutra:\$CUDA_LIB_PATH_64:\$CUDA_LIB_PATH:\$CULA_LIB_PATH_64:\$CULA_LIB_PATH:\$LD_LIBRARY_PATH"
echo
echo "# The following has to be customize according to your Yorick install"
echo "export YORICK_PATH=$YORICK_PATH"
echo "export PATH=\$YORICK_PATH:\$PATH"
echo 'alias yorick="rlwrap -s 2000 yorick"'
