#!/bin/bash

echo "################################################"
echo "ldd -d libcarma"
echo "################################################"
(cd libcarma && ldd -d libcarma.so)
echo "################################################"
#echo "ldd -u libcarma"
#echo "################################################"
#(cd libcarma && ldd -u libcarma.so)
#echo "################################################"
echo "ldd -d libsutra"
echo "################################################"
(cd libsutra && ldd -d libsutra.so)
echo "################################################"
#echo "ldd -u libsutra"
#echo "################################################"
#(cd libsutra && ldd -u libsutra.so)
#echo "################################################"
echo "ldd -d carmaWrap"
echo "################################################"
(cd carmaWrap/build/lib.linux-x86_64-2.7 && ldd -d carmaWrap.so)
echo "################################################"
#echo "ldd -u carmaWrap"
#echo "################################################"
#(cd carmaWrap/build/lib.linux-x86_64-2.7 && ldd -u carmaWrap.so)
#echo "################################################"
echo "ldd -d shesha"
echo "################################################"
(cd shesha/build/lib.linux-x86_64-2.7 && ldd -d shesha.so)
echo "################################################"
#echo "ldd -u shesha"
#echo "################################################"
#(cd shesha/build/lib.linux-x86_64-2.7 && ldd -u shesha.so)
#echo "################################################"
