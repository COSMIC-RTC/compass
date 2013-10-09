(cd libyoga && make clean && make -j) || exit
(make clean install) || exit
(cd yoga_ao/libyoga_ao && make clean && make -j) || exit
(cd yoga_ao && make clean install) || exit
