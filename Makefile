
N_CPU:="$(shell cat /proc/cpuinfo | grep processor | wc -l)"

all: build

.PHONY: build
build:
	python setup.py build_ext

carma: libcarma/libcarma.so

sutra: libsutra/libsutra.so

cython_carmaWrap: carma
	@(cd carmaWrap && make)

cython_shesha: cython_carmaWrap sutra
	@(cd shesha && make)

libcarma/libcarma.so:
	@(cd libcarma && make -j$(N_CPU))

libsutra/libsutra.so: carma
	@(cd libsutra && make -j$(N_CPU))

lib: sutra

cython: cython_shesha

install: all
#	@(cd carmaWrap && make install)
#	@(cd shesha && make install)

clean: clean_carma clean_sutra clean_carmaWrap clean_shesha

clean_carma:
	@(cd libcarma && make clean)

clean_carmaWrap:
	@(cd carmaWrap && make clean)

clean_sutra:
	@(cd libsutra && make clean)

clean_shesha:
	@(cd shesha && make clean)

clean_ao: clean_sutra clean_shesha

clean_lib: clean_carma clean_sutra

clean_cython: clean_carmaWrap clean_shesha

uninstall: clean
	@(cd carmaWrap && make uninstall)
	@(cd shesha && make uninstall)
