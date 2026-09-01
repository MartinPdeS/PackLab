PYTHON ?= python3
BUILD_DIR ?= build
ROOT_DIR := $(CURDIR)
MANUSCRIPT_DIR := docs/manuscript
PYTHON_EXECUTABLE := $(shell $(PYTHON) -c "import sys; print(sys.executable)")
PYBIND11_DIR := $(shell $(PYTHON) -m pybind11 --cmakedir)

.PHONY: configure build install uninstall quick rebuild editable test docs manuscript quality clean

configure:
	cmake -S . -B $(BUILD_DIR) \
		-Dpybind11_DIR="$(PYBIND11_DIR)" \
		-DPython_EXECUTABLE="$(PYTHON_EXECUTABLE)" \
		-DCMAKE_INSTALL_PREFIX="$(ROOT_DIR)"

build:
	cmake --build $(BUILD_DIR) -j

install:
	cmake --install $(BUILD_DIR)
	$(PYTHON) -m pip install --no-deps --no-build-isolation -e .

uninstall:
	$(PYTHON) -m pip uninstall -y PackLab

quick: configure build install

rebuild: configure build install

editable:
	$(PYTHON) -m pip install --no-build-isolation \
		-Cbuild-dir=$(BUILD_DIR) \
		-Ceditable.rebuild=false \
		-Ceditable.mode=inplace \
		-e .

test:
	$(PYTHON) -m pytest

docs:
	$(MAKE) -C docs html SPHINXBUILD="$(PYTHON) -m sphinx"

manuscript:
	MPLBACKEND=Agg $(PYTHON) $(MANUSCRIPT_DIR)/figures/create_figure2.py
	MPLBACKEND=Agg $(PYTHON) $(MANUSCRIPT_DIR)/figures/create_figure3.py
	cd $(MANUSCRIPT_DIR) && latexmk -pdf -interaction=nonstopmode -halt-on-error packlab.tex

quality:
	$(PYTHON) -m flake8 PackLab tests

clean:
	rm -rf $(BUILD_DIR)
