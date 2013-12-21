# these values filled in by    yorick -batch make.i
Y_MAKEDIR=/home/sevin/yorick.git/relocate
Y_EXE=/home/sevin/yorick.git/relocate/bin/yorick
Y_EXE_PKGS=

Y_EXE_HOME=/home/sevin/yorick.git/relocate
Y_EXE_SITE=/home/sevin/yorick.git/relocate
Y_HOME_PKG=

# ----------------------------------------------------- optimization flags

# options for make command line, e.g.-   make COPT=-g TGT=exe
COPT=$(COPT_DEFAULT)
TGT=$(DEFAULT_TGT)
# ------------------------------------------------ macros for this package

PKG_NAME=yoga_template
PKG_I= yoga_template.i

OBJS= yoga_template.o yoga_template.cuo 

# change to give the executable a name other than yorick
#PKG_EXENAME=yorick

MAGMA_DIR=/usr/local/magma
MAGMATEST=/home/sevin/magma-1.4.1
MKLROOT=/opt/intel/composer_xe_2013_sp1/mkl

NVCC = /usr/local/cuda/bin/nvcc 
NVCCFLAGS = -g -DDEBUG -m64 -Xcompiler -fPIC -c --compiler-options -I. -I$(COMPASS_ROOT_DIR)/libcarma/include.h -I$(YOGA_DIR)  -I$(CUDA_INC_PATH) -I$(CULA_INC_PATH) -I$(MAGMA_DIR)/include #-gencode arch=compute_12,code=sm_12  #-I/usr/local/magma-1.4.1/include

# PKG_DEPLIBS=-Lsomedir -lsomelib   for dependencies of this package
PKG_DEPLIBS += $(YOGA_DIR)/yoga.so 
PKG_DEPLIBS += $(MAGMA_DIR)/lib/libmagma.a

# set compiler (or rarely loader) flags specific to this package

PKG_CPPFLAGS = -O3 -DMAGMA_WITH_MKL -DMAGMA_SETAFFINITY -DHAVE_CUBLAS -DADD_ -g -DDEBUG -Wall -fno-common -I$(COMPASS_ROOT_DIR)/libcarma/include.h -I$(YOGA_DIR) -I$(CUDA_INC_PATH) -I$(CULA_INC_PATH) -fPIC -I$(MKLROOT)/include -I$(MAGMA_DIR)/include -I$(MAGMATEST)/control -I$(MAGMATEST)/testing

PKG_LDFLAGS= -L$(MAGMA_DIR)/lib -L$(MKLROOT)/lib/intel64 -L$(CUDA_ROOT)/lib64 -lmagma -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -lcublas -lcudart -lm -lstdc++

# list of additional package names you want in PKG_EXENAME
# (typically Y_EXE_PKGS should be first here)
EXTRA_PKGS=$(Y_EXE_PKGS)

# list of additional files for clean
PKG_CLEAN=

# autoload file for this package, if any
PKG_I_START=
# non-pkg.i include files for this package, if any
I_EXTRA =
PKG_I_EXTRA := $(I_EXTRA)

# -------------------------------- standard targets and rules (in Makepkg)

#--compiler-options -fpermissive 

# set macros Makepkg uses in target and dependency names
# DLL_TARGETS, LIB_TARGETS, EXE_TARGETS
# are any additional targets (defined below) prerequisite to
# the plugin library, archive library, and executable, respectively
PKG_I_DEPS=$(PKG_I)
Y_DISTMAKE=distmake

include $(Y_MAKEDIR)/Make.cfg
include $(Y_MAKEDIR)/Makepkg
include $(Y_MAKEDIR)/Make$(TGT)

# override macros Makepkg sets for rules and other macros
# Y_HOME and Y_SITE in Make.cfg may not be correct (e.g.- relocatable)
Y_HOME=$(Y_EXE_HOME)
Y_SITE=$(Y_EXE_SITE)

# reduce chance of yorick-1.5 corrupting this Makefile
MAKE_TEMPLATE = protect-against-1.5

%.o        : %.cpp
	@printf '\033[36m%s\033[31m%s\033[m\n' "Compiling     " $@
	@$(CC) $(CFLAGS) $(CPPFLAGS) $(PKG_CPPFLAGS) $< -c -o $@

%.cuo      : %.cu
	@printf '\033[36m%s\033[31m%s\033[m\n' "Compiling     " $@
	@$(NVCC) $(NVCCFLAGS) -o $@ $< 

$(PKG_NAME): $(OBJS)
	@printf '\033[36m%s\033[31m%s\033[m\n' "Linking       " $(PKG_NAME)
	@$(CC) $(PKG_CPPFLAGS) $(PKG_DEPLIBS) -o $@ $^ 

#configure:
#	yorick -batch make.i

# -------------------------------------------------------- end of Makefile
