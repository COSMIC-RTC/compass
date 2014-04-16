
N_CPU:=$(cat /proc/cpuinfo | grep processor | wc -l)

all: configure
	(cd libcarma && make -j$(N_CPU))
	(cd yoga && make install)
	(cd libsutra && make -j$(N_CPU))
	(cd yoga_ao && make install)

clean: configure
	(cd libcarma && make clean)
	(cd yoga && make clean)
	(cd libsutra && make clean)
	(cd yoga_ao && make clean)

configure:
	./configure.sh

