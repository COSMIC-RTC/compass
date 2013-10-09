(cd libyoga && make -j) || exit
(make clean install) || exit
(cd yoga_ao/libyoga_ao && make -j) || exit
(cd yoga_ao && make clean install) || exit
