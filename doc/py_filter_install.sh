git clone https://github.com/Feneric/doxypypy

(cd doxypypy; python -m setup install)

cp $COMPASS_ROOT/doc/py_filter $(dirname $(which python))
