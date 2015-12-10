all: install

build:
	python setup.py build_ext

clean:
	python clean.py && rm -rf build

install: build
	python setup.py install --user

uninstall:
	$(error "not handle yet")
