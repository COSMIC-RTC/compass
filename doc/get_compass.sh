#!/bin/bash

# wget -O- https://lesia.obspm.fr/compass/get_compass.sh | /bin/bash
# curl -L https://lesia.obspm.fr/compass/get_compass.sh | /bin/bash

check_linux() {
if [ -f /etc/os-release ]; then
    # freedesktop.org and systemd
    . /etc/os-release
    OS=$NAME
    VER=$VERSION_ID
elif type lsb_release >/dev/null 2>&1; then
    # linuxbase.org
    OS=$(lsb_release -si)
    VER=$(lsb_release -sr)
elif [ -f /etc/lsb-release ]; then
    # For some versions of Debian/Ubuntu without lsb_release command
    . /etc/lsb-release
    OS=$DISTRIB_ID
    VER=$DISTRIB_RELEASE
elif [ -f /etc/debian_version ]; then
    # Older Debian/Ubuntu/etc.
    OS=Debian
    VER=$(cat /etc/debian_version)
elif [ -f /etc/SuSe-release ]; then
    # Older SuSE/etc.
    ...
elif [ -f /etc/redhat-release ]; then
    # Older Red Hat, CentOS, etc.
    ...
else
    # Fall back to uname, e.g. "Linux <version>", also works for BSD, etc.
    OS=$(uname -s)
    VER=$(uname -r)
fi

	echo "linux detected: $OS $VER"
	if [ "x$OS" == "xCentOS Linux" ] && [ "x$VER" == "x7" ]
	then
		echo "special binaries for CentOS 7"
		EXTENSION=".centos7"
	fi
}

check_cuda_version() {
    nvcc --version > /dev/null

    if [ $? != 0 ]
    then
        echo 'Please install CUDA Toolkit'
        debug_conf
    fi

   CUDA_VERSION=`nvcc --version | grep release | awk -F, '{ print $2 }' | awk '{print $2}' | sed 's:\.::'`
}

check_conda() {
    conda --version > /dev/null
    if [ $? != 0 ]
    then
        cat >> $HOME/.bashrc << OEF
## CONDA default definitions
export CONDA_ROOT=$HOME/miniconda3
export PATH=$CONDA_ROOT/bin:$PATH
OEF
		export CONDA_ROOT=$HOME/miniconda3
		export PATH=$CONDA_ROOT/bin:$PATH
        mkdir -p $HOME/tmp_compass
        cd $HOME/tmp_compass
        wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
        bash Miniconda3-latest-Linux-x86_64.sh -b -p $CONDA_ROOT
    else
        echo "conda already installed, skipping"
    fi
}

get_compass() {
	export SHESHA_ROOT=$HOME/shesha
	export PYTHONPATH=$NAGA_ROOT:$SHESHA_ROOT:$PYTHONPATH
    cat >> $HOME/.bashrc << OEF
#COMPASS default definitions
export SHESHA_ROOT=$HOME/shesha
export PYTHONPATH=$NAGA_ROOT:$SHESHA_ROOT:$PYTHONPATH
OEF
	export SHESHA_ROOT=$HOME/shesha
	export PYTHONPATH=$NAGA_ROOT:$SHESHA_ROOT:$PYTHONPATH
	echo "install compass from label $CUDA_VERSION$EXTENSION"
    conda install -y -c compass/label/cuda$CUDA_VERSION$EXTENSION compass aenum
    cd $HOME
    git clone https://github.com/ANR-COMPASS/shesha.git
}

test_compass() {
    cd $SHESHA_ROOT
    ipython -i shesha/scripts/closed_loop.py data/par/par4bench/scao_sh_16x16_8pix.py
}

debug_conf () {
	echo " ############################################################"
	echo " #### DEBUG information to send at arnaud.sevin@obspm.fr ####"
	echo " ############################################################"
	echo
	check_linux
	echo
	cat /etc/*release
	echo
	echo "lspci | grep -i vga"
	lspci | grep -i vga
	echo
	echo "lspci | grep -i nvidia"
	lspci | grep -i nvidia
	echo
	echo "nvcc version"
	nvcc --version
	echo
	echo "nvcc path"
	which nvcc
	echo
	echo "list $HOME"
    ls -C1 $HOME
	echo
	exit 1
}

check_linux
check_cuda_version
check_conda || debug_conf
get_compass || debug_conf
source $HOME/.bashrc
test_compass || debug_conf
