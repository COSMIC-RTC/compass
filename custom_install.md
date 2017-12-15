- [COMPASS installation: the hard way](#compass-installation-the-hard-way)
    - [OpenBLAS installation](#openblas-installation)
    - [MAGMA installation](#magma-installation)
    - [CONDA minimal installation](#conda-minimal-installation)

# COMPASS installation: the hard way

## OpenBLAS installation

```bash
git clone https://github.com/xianyi/OpenBLAS.git
cd OpenBLAS/
make USE_OPENMP=1
make install PREFIX=$HOME/local/openblas
export OPENBLAS_ROOT=$HOME/local/openblas
```

## MAGMA installation

```bash
wget http://icl.utk.edu/projectsfiles/magma/downloads/magma-2.3.0.tar.gz
tar xf magma-2.3.0.tar.gz
cd magma-2.3.0
cp make.inc-examples/make.inc.openblas make.inc
GPU_TARGET=sm_52 OPENBLASDIR=$OPENBLAS_ROOT CUDADIR=/opt/cuda make -j 8 shared sparse-shared
GPU_TARGET=sm_52 OPENBLASDIR=$OPENBLAS_ROOT CUDADIR=/opt/cuda make install prefix=$HOME/local/magma
```

## CONDA minimal installation

```bash
export CONDA_ROOT=$HOME/miniconda3
bash Miniconda3-latest-Linux-x86_64.sh -b -p $CONDA_ROOT  

pip install cython

git clone https://github.com/numpy/numpy.git numpy
cd numpy/
cp site.cfg.example site.cfg
```

[fill openblas part in site.cgf]

```bash
python setup.py build -j 8 install --prefix $CONDA_ROOT
cd  $CONDA_ROOT/lib/python3.6/site-packages
unzip numpy-1.14.0.dev0+0058e70-py3.6-linux-x86_64.egg

pip install pyqt5 ipython astropy blaze h5py nose scipy pyqtgraph matplotlib tqdm docopt
```
