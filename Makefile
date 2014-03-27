all: configure
	(cd libcarma && make -j12)
	(cd yoga && make install)
	(cd libsutra && make -j12)
	(cd yoga_ao && make install)

clean: configure
	(cd libcarma && make clean)
	(cd yoga && make clean)
	(cd libsutra && make clean)
	(cd yoga_ao && make clean)

configure:
	./configure.sh

