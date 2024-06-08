.PHONY: dist examples

all: examples

test: 
	python -m pytest

install-dev:
	python -m pip install -e .

install:
	python -m pip install .

examples:
	cd examples && make

cleanall:
	cd examples && make cleanall

dist:
	python -m build
