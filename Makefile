.PHONY: dist examples

all: test

test: 
	python -m pytest

install-dev:
	python -m pip install -e .

install:
	python -m pip install .

examples:
	cd examples && make

cleanrun:
	cd examples && make cleanall

dist:
	python -m build
