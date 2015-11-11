
N_CPU:="$(shell cat /proc/cpuinfo | grep processor | wc -l)"

YORICK = $(shell which yorick)
YORICK-URL="https://github.com/dhmunro/yorick.git"

check_yorick: 
ifeq ($(YORICK),)
	@echo "yorick not in PATH"
	@echo "please install yorick with the commmand : make yorick"
	@echo "and add the relocate/bin directory in the PATH"
	@false
endif

yorick:
	@printf '\033[36m%s\033[31m%s\033[m\n' "Compiling     " $@
	@git clone $(YORICK-URL) yorick.git
	@(cd yorick.git && ./configure && make -j && make install)
	@echo "REMINDER : put this line into your .bashrc"
	@echo "export PATH=`pwd`/yorick.git/relocate/bin:\$$PATH"

all: configure
	@echo "using $(N_CPU) jobs"
	@(if [ -z "$(COMPILATION_LAM)" ]; then echo "lam won't be compiled." ; else cd lam && make; fi)
	@(cd libcarma && make -j$(N_CPU))
	@(cd yoga && make)
	@(cd libsutra && make -j$(N_CPU))
	@(cd yoga_ao && make)

install: all configure
	@(cd yoga && make install)
	@(cd yoga_ao && make install)

clean: configure
	@(cd libcarma && make clean)
	@(cd yoga && make clean)
	@(cd libsutra && make clean)
	@(cd yoga_ao && make clean)
	@(if [ -z "$(COMPILATION_LAM)" ]; then echo "lam has been ignored." ; else cd lam && make clean; fi)

uninstall: configure clean
	@(cd yoga && make uninstall)
	@(cd yoga_ao && make uninstall)

configure: check_yorick
	./configure.sh

