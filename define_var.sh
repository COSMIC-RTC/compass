# uncomment after editing the file
echo "Please edit this file before run it"

# CUDA default definitions
#export GENCODE='"arch=compute_12,code=sm_12" -DGPUSHMEM=120' #for old tesla
#export GENCODE='"arch=compute_20,code=sm_20" -DGPUSHMEM=200' #for Fermi
#export GENCODE='"arch=compute_30,code=sm_30" -DGPUSHMEM=300' #for Kepler
#export GENCODE='"arch=compute_35,code=sm_35" -DGPUSHMEM=300' #for K40
export GENCODE='"arch=compute_37,code=sm_37" -DGPUSHMEM=300' #for K40
#${GENCODE?"Need to set GENCODE"}

export CUDA_ROOT= #/usr/local/cuda
export CUDA_INC_PATH=$CUDA_ROOT/include
export CUDA_LIB_PATH=$CUDA_ROOT/lib
export CUDA_LIB_PATH_64=$CUDA_ROOT/lib64
export CPLUS_INCLUDE_PATH=$CUDA_INC_PATH
export PATH=$CUDA_ROOT/bin:$PATH

# CULA default definitions
# export CULA_ROOT= #/usr/local/cula
# export CULA_INC_PATH= #$CULA_ROOT/include
# export CULA_LIB_PATH= #$CULA_ROOT/lib
# export CULA_LIB_PATH_64= #$CULA_ROOT/lib64

# MAGMA definitions (uncomment this line if MAGMA is installed)
export MAGMA_ROOT= #$HOME/local/cuda
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MAGMA_ROOT/lib
export PKG_CONFIG_PATH=$MAGMA_ROOT/lib/pkgconfig

# COMPASS default definitions
export COMPASS_ROOT= #$HOME/compass
export NAGA_ROOT=$COMPASS_ROOT/carmaWrap
export SHESHA_ROOT=$COMPASS_ROOT/shesha

export LD_LIBRARY_PATH=$COMPASS_ROOT/libcarma:$COMPASS_ROOT/libsutra:$CUDA_LIB_PATH_64:$CUDA_LIB_PATH:$CULA_LIB_PATH_64:$CULA_LIB_PATH:$LD_LIBRARY_PATH

export PYTHONPATH=$SHESHA/src:$NAGA/src:$SHESHA/lib:$NAGA/lib:$PYTHONPATH

echo "COMPASS will use this configuration, please add those lines into the .bashrc"

echo "# CUDA default definitions"
echo "# The value of GENCODE should be modified according to your GPU generation"
echo "export CUDA_ROOT=$CUDA_ROOT"
echo "export CUDA_INC_PATH=\$CUDA_ROOT/include"
echo "export CUDA_LIB_PATH=\$CUDA_ROOT/lib"
echo "export CUDA_LIB_PATH_64=\$CUDA_ROOT/lib64"
echo "export CPLUS_INCLUDE_PATH=\$CUDA_INC_PATH"
echo "export PATH=\$CUDA_ROOT/bin:\$PATH"
echo
echo "export GENCODE=\"$GENCODE\""
echo
# echo "# CULA default definitions (optional)"
# echo "export CULA_ROOT=$CULA_ROOT"
# echo "export CULA_INC_PATH=$CULA_INC_PATH"
# echo "export CULA_LIB_PATH=$CULA_LIB_PATH"
# echo "export CULA_LIB_PATH_64=$CULA_LIB_PATH_64"
# echo
echo "# MAGMA definitions (optional)"
echo "export MAGMA_ROOT=$MAGMA_ROOT"
echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:\$MAGMA_ROOT/lib"
echo "export PKG_CONFIG_PATH=\$PKG_CONFIG_PATH:\$MAGMA_ROOT/lib/pkgconfig"
echo
echo "# COMPASS default definitions"
echo "export COMPASS_ROOT=$COMPASS_ROOT"
echo "export NAGA_ROOT=\$COMPASS_ROOT/carmaWrap"
echo "export SHESHA_ROOT=\$COMPASS_ROOT/shesha"
echo
echo "export LD_LIBRARY_PATH=\$COMPASS_ROOT/libcarma:\$COMPASS_ROOT/libsutra:\$CUDA_LIB_PATH_64:\$CUDA_LIB_PATH:\$LD_LIBRARY_PATH"
echo "export PYTHONPATH=\$SHESHA_ROOT/src:\$NAGA_ROOT/src:\$SHESHA_ROOT/lib:\$NAGA_ROOT/lib:\$PYTHONPATH"
