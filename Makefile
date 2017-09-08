
N_CPU:="$(shell cat /proc/cpuinfo | grep processor | wc -l)"

all:
	@echo "using $(N_CPU) jobs"
	@(if [ -z "$(COMPILATION_LAM)" ]; then echo "lam won't be compiled." ; else cd lam && make; fi)
	@(cd libcarma && make -j$(N_CPU))
	@(cd naga && make)
	@(cd libsutra && make -j$(N_CPU))
	@(cd shesha && make)

install: all
#	@(cd naga && make install)
#	@(cd shesha && make install)

clean:
	@(cd libcarma && make clean)
	@(cd naga && make clean)
	@(cd libsutra && make clean)
	@(cd shesha && make clean)
	@(if [ -z "$(COMPILATION_LAM)" ]; then echo "lam has been ignored." ; else cd lam && make clean; fi)

clean_ao:
	@(cd libsutra && make clean)
	@(cd shesha && make clean)
	@(if [ -z "$(COMPILATION_LAM)" ]; then echo "lam has been ignored." ; else cd lam && make clean; fi)

clean_lib:
	@(cd libsutra && make clean)
	@(cd libcarma && make clean)

clean_cython:
	@(cd naga && make clean)
	@(cd shesha && make clean)

lib:
	@(cd libcarma && make -j$(N_CPU))
	@(cd libsutra && make -j$(N_CPU))

cython:
	@(cd naga && make)
	@(cd shesha && make)

uninstall: clean
	@(cd naga && make uninstall)
	@(cd shesha && make uninstall)
