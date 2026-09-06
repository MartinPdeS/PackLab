.DEFAULT_GOAL := help

PYTHON ?= python3
BUILD_DIR ?= build
ROOT_DIR := $(CURDIR)
MANUSCRIPT_DIR ?= docs/manuscript/softwareX
MANUSCRIPT_TEX ?= packlab.tex
REGENERATE_FIGURES ?= 0
TAG_VERSION ?= $(or $(VERSION),$(filter v%,$(MAKECMDGOALS)))
RELEASE_KIND := $(filter major minor patch,$(MAKECMDGOALS))
PYTHON_EXECUTABLE = $(shell $(PYTHON) -c "import sys; print(sys.executable)")
PYBIND11_DIR = $(shell $(PYTHON) -m pybind11 --cmakedir)
MPLBACKEND ?= Agg
MPLCONFIGDIR ?= $(abspath $(BUILD_DIR))/matplotlib
XDG_CACHE_HOME ?= $(abspath $(BUILD_DIR))/cache
export MPLBACKEND MPLCONFIGDIR XDG_CACHE_HOME
MANUSCRIPT_FIGURES := \
	$(MANUSCRIPT_DIR)/figures/create_figure2.py \
	$(MANUSCRIPT_DIR)/figures/create_figure3.py \
	$(MANUSCRIPT_DIR)/figures/create_figure4_independent_validation.py
LINT_PATHS := PackLab tests tools development/compare_mctpy.py \
	development/extract_mctpy_reference.py $(MANUSCRIPT_DIR)/figures
FORMAT_PATHS := tools/check_manuscript.py tools/check_release.py tools/clean.py \
	tools/project_doctor.py tests/test_workflow_tools.py \
	development/compare_mctpy.py development/extract_mctpy_reference.py \
	$(MANUSCRIPT_DIR)/figures

.PHONY: help setup workflow-dirs configure build install uninstall quick rebuild editable \
	test test-fast docs manuscript manuscript-figures manuscript-validation \
	manuscript-check manuscript-clean reproduce-paper quality check doctor \
	release-check tag release major minor patch clean

.NOTPARALLEL: check reproduce-paper

help:
	@echo "PackLab development commands"
	@echo ""
	@echo "  make setup                 Install editable development dependencies"
	@echo "  make build                 Build the native extensions"
	@echo "  make test                  Run the full test suite with coverage"
	@echo "  make test-fast             Run tests without coverage reporting"
	@echo "  make quality               Check formatting and linting"
	@echo "  make docs                  Build the Sphinx documentation"
	@echo "  make manuscript            Build the SoftwareX manuscript"
	@echo "  make manuscript-figures    Regenerate quantitative manuscript figures"
	@echo "  make manuscript-validation Compare PackLab with stored reference data"
	@echo "  make manuscript-check      Build and check SoftwareX submission limits"
	@echo "  make reproduce-paper       Recreate validation, figures, and manuscript"
	@echo "  make check                 Run local pre-push checks"
	@echo "  make doctor                Diagnose the local development environment"
	@echo "  make release-check         Check version metadata consistency"
	@echo "  make tag VERSION=vX.Y.Z    Create a release commit and annotated tag"
	@echo "  make release patch         Create the next patch release and tag"
	@echo "  make release minor         Create the next minor release and tag"
	@echo "  make release major         Create the next major release and tag"
	@echo "  make clean                 Remove generated build and manuscript files"

setup:
	$(PYTHON) -m pip install -e ".[testing,documentation,dev,scattering]"

workflow-dirs:
	$(PYTHON) -c "from pathlib import Path; [Path(path).mkdir(parents=True, exist_ok=True) for path in ('$(MPLCONFIGDIR)', '$(XDG_CACHE_HOME)')]"

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

editable: configure build install

test: workflow-dirs
	$(PYTHON) -m pytest

test-fast: workflow-dirs
	$(PYTHON) -m pytest -o addopts="-q"

docs:
	$(MAKE) -C docs html SPHINXBUILD="$(PYTHON_EXECUTABLE) -m sphinx"

manuscript:
ifeq ($(REGENERATE_FIGURES),1)
	$(MAKE) manuscript-figures
endif
	cd $(MANUSCRIPT_DIR) && latexmk -g -pdf -interaction=nonstopmode -halt-on-error $(MANUSCRIPT_TEX)

manuscript-figures: workflow-dirs
	$(PYTHON) $(word 1,$(MANUSCRIPT_FIGURES))
	$(PYTHON) $(word 2,$(MANUSCRIPT_FIGURES))
	$(PYTHON) $(word 3,$(MANUSCRIPT_FIGURES))

manuscript-validation: workflow-dirs
	$(PYTHON) -c "from pathlib import Path; Path('$(BUILD_DIR)/manuscript-validation').mkdir(parents=True, exist_ok=True)"
	$(PYTHON) \
		development/compare_mctpy.py \
		--output-directory $(BUILD_DIR)/manuscript-validation

manuscript-check: manuscript
	$(PYTHON) tools/check_manuscript.py --manuscript-dir $(MANUSCRIPT_DIR)

manuscript-clean:
	$(PYTHON) tools/clean.py --manuscript-only --manuscript-dir $(MANUSCRIPT_DIR)

reproduce-paper:
	$(MAKE) manuscript-figures
	$(MAKE) manuscript-validation
	$(MAKE) manuscript-check

quality:
	$(PYTHON) -m ruff check --ignore E501 $(LINT_PATHS)
	$(PYTHON) -m ruff format --check $(FORMAT_PATHS)

check:
	$(MAKE) editable
	$(MAKE) quality
	$(MAKE) test
	$(MAKE) docs
	$(MAKE) manuscript-check

doctor:
	$(PYTHON) tools/project_doctor.py

release-check:
	$(PYTHON) tools/check_release.py $(if $(VERSION),--version $(VERSION),)

# Create a release commit and annotated tag. Both forms are supported:
#   make tag v0.7.1
#   make tag VERSION=v0.7.1
tag:
	$(PYTHON) tools/release_tag.py "$(TAG_VERSION)"

# Derive the next tag from the highest existing vMAJOR.MINOR.PATCH release.
# Examples: make release patch, make release minor, make release major
release:
	@test "$(words $(RELEASE_KIND))" -eq 1 || { echo "usage: make release [patch|minor|major]" >&2; exit 2; }
	$(PYTHON) tools/release_tag.py "$$($(PYTHON) tools/next_release_version.py $(RELEASE_KIND))"

major minor patch:
	@:

clean:
	$(PYTHON) tools/clean.py --build-dir "$(BUILD_DIR)" --build-dir .skbuild \
		--manuscript-dir $(MANUSCRIPT_DIR)
