#!/bin/env bash

echo -e "# COMPASS configuration script output"
echo -e "\n## hostname\n"
hostname
echo -e "\n## lspci\n"
lspci | grep -i vga
echo -e "\n## OS\n"
cat /etc/lsb-release
echo -e "\n## kernel\n"
uname -a
echo -e "\n## shesha version\n"
ipython -c 'import shesha; print(shesha.__version__)'
echo -e "\n## conda list\n"
conda list

CONDA=`conda list | grep compass | wc -l`

if [[ $CONDA -eq "0" ]]
then
	echo -e "\n## development version\n"
	echo "ldd -d libcarma"
        ldd -d local/lib/libcarma.so
	echo "ldd -d libsutra"
        ldd -d local/lib/libsutra.so
	echo "ldd -d carmaWrap"
        ldd -d local/python/carmaWrap.cpython-*.so | grep -v Py
	echo "ldd -d sutrarap"
        ldd -d local/python/sutraWrap.cpython-*.so  | grep -v Py
else
	echo -e "\n## public version\n"
	conda list | grep compass
fi
