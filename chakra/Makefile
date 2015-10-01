all: install

build:
	python setup.py build_ext -i

clean:
	python clean.py

install: build
	python setup.py install --user
