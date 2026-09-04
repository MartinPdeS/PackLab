"""Tests for repository workflow helpers."""

import importlib.util
from pathlib import Path
from types import ModuleType

import pytest


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


def load_tool_module(name: str) -> ModuleType:
    """Load a repository-only workflow helper without packaging it."""
    path = REPOSITORY_ROOT / "tools" / f"{name}.py"
    specification = importlib.util.spec_from_file_location(f"packlab_test_{name}", path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"could not load workflow helper: {path}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


check_manuscript = load_tool_module("check_manuscript")
check_release = load_tool_module("check_release")
clean = load_tool_module("clean")


def test_extract_latex_environment() -> None:
    """Manuscript checks extract the requested environment only."""
    source = r"\begin{keyword}one \sep two\end{keyword}"
    assert check_manuscript.extract_environment(source, "keyword") == r"one \sep two"


def test_release_versions_are_aligned() -> None:
    """Package release artifacts declare one version."""
    declared = check_release.versions()
    assert len(set(declared.values())) == 1


@pytest.mark.parametrize(
    ("version", "expected"),
    [("1.2.3", "1.2.3"), ("v1.2.3", "1.2.3")],
)
def test_normalized_version(version: str, expected: str) -> None:
    """Release checks accept plain and Git-style version strings."""
    assert check_release.normalized_version(version) == expected


def test_cleaner_rejects_repository_root() -> None:
    """The cleaning helper cannot recursively remove the repository."""
    with pytest.raises(ValueError, match="outside the repository"):
        clean.safe_directory(clean.ROOT)


def test_cleaner_accepts_scoped_directory() -> None:
    """A generated directory inside the repository is a valid clean target."""
    build_directory = clean.ROOT / "build" / "workflow-test"
    assert clean.safe_directory(build_directory) == Path(build_directory).resolve()
