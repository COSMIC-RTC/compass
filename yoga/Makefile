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

DEBUG=1

I_DIR=yorick

PKG_NAME = yoga
PKG_I=$(I_DIR)/yoga.i

OBJS=yoga.o 

#MAGMA_PATH= /home/sevin/magma-1.4.1

# change to give the executable a name other than yorick
#PKG_EXENAME=yorick

# PKG_DEPLIBS=-Lsomedir -lsomelib   for dependencies of this package
PKG_DEPLIBS=  -L$(YOGA_DIR) -L$(COMPASS_ROOT_DIR)/libcarma -lcarma -lstdc++ # 

# set compiler (or rarely loader) flags specific to this package
#PKG_CPPFLAGS= -Wall -g -O3 -DUNIX -DADD_ -DGPUSHMEM=130 -fPIC -Xlinker -zmuldefs -fno-common -I$(H_DIR) -I/usr/local/cuda/include 
PKG_CPPFLAGS= -Wall -fno-common -I$(COMPASS_ROOT_DIR)/libcarma/include.h -Iinclude.h -I$(CUDA_INC_PATH) -I$(CULA_INC_PATH) -DUSE_MAGMA
ifneq ($(MAGMA_PATH),)
 ifeq ($(MAGMA_USE),ATLAS)
  PKG_CPPFLAGS +=-I$(MAGMA_PATH)/include -DUSE_MAGMA -DHAVE_CUBLAS -DADD_ -DMAGMA_SETAFFINITY
 endif
endif

ifneq ($(DEBUG),)
PKG_CPPFLAGS+= -g -DDEBUG
endif

#-I/usr/local/cula/include
PKG_LDFLAGS=

# list of additional package names you want in PKG_EXENAME
# (typically Y_EXE_PKGS should be first here)
EXTRA_PKGS=$(Y_EXE_PKGS)

# list of additional files for clean
PKG_CLEAN=

# autoload file for this package, if any
PKG_I_START=
# non-pkg.i include files for this package, if any
PKG_I_EXTRA= $(I_DIR)/check_yoga.i $(I_DIR)/fits-utils.i $(I_DIR)/util_fr.i
# yoga_fft.i yoga_matmult.i yoga_utils.i

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

# ------------------------------------- targets and rules for this package
# color for printf
#31=red, 32=green, 33=yellow,34=blue, 35=pink, 36=cyan, 37=white

yoga.o: yoga.cpp
	@printf '\033[36m%s\033[31m%s\033[m\n' "Compiling     " $@
	$(CC) $(CFLAGS) $(PKG_CPPFLAGS) $< -c -o $@

$(PKG_NAME): $(OBJS)
	@printf '\033[36m%s\033[31m%s\033[m\n' "Linking       " $(PKG_NAME)
	$(CC) $(PKG_DEPLIBS) -o $@ $^ 

