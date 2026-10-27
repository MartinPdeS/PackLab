"""Tests for repository workflow helpers."""

import importlib.util
from importlib.metadata import version as distribution_version
from pathlib import Path
import runpy
from types import ModuleType

from packaging.version import Version
import pytest

from martinpdes_maintenance.cli import load, versions


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
TOOLS_DIRECTORY = REPOSITORY_ROOT / "tools"

if not TOOLS_DIRECTORY.is_dir():
    pytest.skip(
        "repository workflow helpers are not staged in this package test environment",
        allow_module_level=True,
    )


def load_tool_module(name: str) -> ModuleType:
    """Load a repository-only workflow helper without packaging it."""
    path = TOOLS_DIRECTORY / f"{name}.py"
    specification = importlib.util.spec_from_file_location(f"packlab_test_{name}", path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"could not load workflow helper: {path}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


check_manuscript = load_tool_module("check_manuscript")
clean = load_tool_module("clean")


def test_extract_latex_environment() -> None:
    """Manuscript checks extract the requested environment only."""
    source = r"\begin{keyword}one \sep two\end{keyword}"
    assert check_manuscript.extract_environment(source, "keyword") == r"one \sep two"


def test_release_versions_are_aligned() -> None:
    """Package release artifacts declare one version."""
    _, config = load(REPOSITORY_ROOT)
    declared = versions(REPOSITORY_ROOT, config)
    assert len(set(declared.values())) == 1


def test_docs_version_uses_package_version(monkeypatch) -> None:
    """The Sphinx version and switcher use _version.py, even when a tag is set."""
    pytest.importorskip("martinpdes_docs")
    if Version(distribution_version("MartinPdeS-docs")) < Version("0.4.0"):
        pytest.skip("automatic version discovery requires MartinPdeS-docs 0.4.0")
    monkeypatch.setenv("tag", "v9.9.9")
    config = runpy.run_path(str(REPOSITORY_ROOT / "docs" / "source" / "conf.py"))
    expected = runpy.run_path(str(REPOSITORY_ROOT / "PackLab" / "_version.py"))["__version__"]
    assert config["version"] == config["release"] == expected
    assert config["html_theme_options"]["switcher"]["version_match"] == expected


def test_cleaner_rejects_repository_root() -> None:
    """The cleaning helper cannot recursively remove the repository."""
    with pytest.raises(ValueError, match="outside the repository"):
        clean.safe_directory(clean.ROOT)


def test_cleaner_accepts_scoped_directory() -> None:
    """A generated directory inside the repository is a valid clean target."""
    build_directory = clean.ROOT / "build" / "workflow-test"
    assert clean.safe_directory(build_directory) == Path(build_directory).resolve()
