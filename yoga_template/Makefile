# these values filled in by    yorick -batch make.i
Y_MAKEDIR=/home/brujo/yorick-2.2/relocate
Y_EXE=/home/brujo/yorick-2.2/relocate/bin/yorick
Y_EXE_PKGS=

Y_EXE_HOME=/home/brujo/yorick-2.2/relocate
Y_EXE_SITE=/home/brujo/yorick-2.2/relocate
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

NVCC = /usr/local/cuda/bin/nvcc 
NVCCFLAGS = -g -DDEBUG -m64 -Xcompiler -fPIC -c --compiler-options -I. -I$(YOGA_DIR)/libcarma/include.h -I$(YOGA_DIR) -I/usr/local/cuda/include -I/usr/local/cula/include -gencode arch=compute_12,code=sm_12 

# PKG_DEPLIBS=-Lsomedir -lsomelib   for dependencies of this package
PKG_DEPLIBS += $(YOGA_DIR)/yoga.so

# set compiler (or rarely loader) flags specific to this package

PKG_CPPFLAGS = -g -DDEBUG -Wall -fno-common -I$(YOGA_DIR)/libcarma/include.h -I$(YOGA_DIR) -I/usr/local/cuda/include -I/usr/local/cula/include -fPIC

PKG_LDFLAGS=  -lstdc++

# list of additional package names you want in PKG_EXENAME
# (typically Y_EXE_PKGS should be first here)
EXTRA_PKGS=$(Y_EXE_PKGS)

# list of additional files for clean
PKG_CLEAN=

# autoload file for this package, if any
PKG_I_START=
# non-pkg.i include files for this package, if any
I_EXTRA = rico.i
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
