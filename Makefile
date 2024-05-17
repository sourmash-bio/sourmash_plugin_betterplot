all: test

test: 
	python -m pytest

install-dev:
	python -m pip install -e .

install:
	python -m pip install .
