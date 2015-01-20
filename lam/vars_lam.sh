#DIR=`realpath $0`
#DIR=`dirname $DIR`
#export COMPASS_ROOT_DIR=`dirname $DIR`
export COMPASS_ROOT_DIR=/home/tgautrais/compass_var_lam/trunk
export yorick="$COMPASS_ROOT_DIR/yoga/yorick"
export YOGA_DIR="$COMPASS_ROOT_DIR/yoga"
export YOGA_AO_DIR="$COMPASS_ROOT_DIR/yoga_ao"
export LIBCARMA_DIR="$COMPASS_ROOT_DIR/libcarma"
export LIBSUTRA_DIR="$COMPASS_ROOT_DIR/libsutra"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$YOGA_DIR
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LIBCARMA_DIR
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LIBSUTRA_DIR
export YOGA_AO_TOP=$YOGA_AO_DIR

export COMPILATION_LAM='  '
