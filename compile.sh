(cd libcarma && make -j8) || exit
(cd yoga && make clean install) || exit
(cd libsutra && make -j8) || exit
(cd yoga_ao && make clean install) || exit
