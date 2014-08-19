
N_CPU:="$(shell cat /proc/cpuinfo | grep processor | wc -l)"

all: configure
	@echo "using $(N_CPU) jobs"
	@(if [ -z "$(COMPILATION_LAM)" ]; then echo "lam/kalman_CPU_GPU won't be compiled." ; else cd lam/kalman_CPU_GPU && make; fi)
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
	@(if [ -z "$(COMPILATION_LAM)" ]; then echo "lam/kalman_CPU_GPU has been ignored." ; else cd lam/kalman_CPU_GPU && make clean; fi)



uninstall: configure clean
	@(cd yoga && make uninstall)
	@(cd yoga_ao && make uninstall)

configure:
	./configure.sh

