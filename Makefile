PYTHON ?= python3
BUILD_DIR ?= build
ROOT_DIR := $(CURDIR)
MANUSCRIPT_DIR := docs/manuscript
REGENERATE_FIGURES ?= 0
TAG_VERSION ?= $(or $(VERSION),$(filter v%,$(MAKECMDGOALS)))
PYTHON_EXECUTABLE = $(shell $(PYTHON) -c "import sys; print(sys.executable)")
PYBIND11_DIR = $(shell $(PYTHON) -m pybind11 --cmakedir)

.PHONY: configure build install uninstall quick rebuild editable test docs manuscript quality tag clean

ifneq ($(filter tag,$(MAKECMDGOALS)),)
ifneq ($(strip $(TAG_VERSION)),)
.PHONY: $(TAG_VERSION)
$(TAG_VERSION):
	@:
endif
endif

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

ifeq ($(REGENERATE_FIGURES),1)
	MPLBACKEND=Agg $(PYTHON) $(MANUSCRIPT_DIR)/figures/create_figure2.py
	MPLBACKEND=Agg $(PYTHON) $(MANUSCRIPT_DIR)/figures/create_figure3.py
	MPLBACKEND=Agg $(PYTHON) $(MANUSCRIPT_DIR)/figures/create_figure4_independent_validation.py
endif
	cd $(MANUSCRIPT_DIR) && latexmk -pdf -interaction=nonstopmode -halt-on-error packlab.tex

quality:
	$(PYTHON) -m ruff check PackLab tests
	$(PYTHON) -m ruff format --check PackLab tests

# Create a release commit and annotated tag. Both forms are supported:
#   make tag v0.5.0
#   make tag VERSION=v0.5.0
tag:
	$(PYTHON) tools/release_tag.py "$(TAG_VERSION)"

clean:
	rm -rf $(BUILD_DIR)
