(cd libcarma && make -j) || exit
(make clean install) || exit
(cd libsutra && make -j) || exit
(cd yoga_ao && make clean install) || exit
