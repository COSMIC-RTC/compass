
N_CPU:="$(shell cat /proc/cpuinfo | grep processor | wc -l)"

all: build

.PHONY: build
build:
	python setup.py build_ext

carma: libcarma/libcarma.so

sutra: libsutra/libsutra.so

cython_naga: carma
	@(cd naga && make)

cython_shesha: cython_naga sutra
	@(cd shesha && make)

libcarma/libcarma.so:
	@(cd libcarma && make -j$(N_CPU))

libsutra/libsutra.so: carma
	@(cd libsutra && make -j$(N_CPU))

lib: sutra

cython: cython_shesha

install: all
#	@(cd naga && make install)
#	@(cd shesha && make install)

clean: clean_carma clean_sutra clean_naga clean_shesha

clean_carma:
	@(cd libcarma && make clean)

clean_naga:
	@(cd naga && make clean)

clean_sutra:
	@(cd libsutra && make clean)

clean_shesha:
	@(cd shesha && make clean)

clean_ao: clean_sutra clean_shesha

clean_lib: clean_carma clean_sutra

clean_cython: clean_naga clean_shesha

uninstall: clean
	@(cd naga && make uninstall)
	@(cd shesha && make uninstall)
