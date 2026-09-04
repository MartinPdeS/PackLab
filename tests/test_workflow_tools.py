"""Tests for repository workflow helpers."""

from pathlib import Path

import pytest

from tools.check_manuscript import extract_environment
from tools.check_release import normalized_version, versions
from tools.clean import ROOT, safe_directory


def test_extract_latex_environment() -> None:
    """Manuscript checks extract the requested environment only."""
    source = r"\begin{keyword}one \sep two\end{keyword}"
    assert extract_environment(source, "keyword") == r"one \sep two"


def test_release_versions_are_aligned() -> None:
    """Package release artifacts declare one version."""
    declared = versions()
    assert len(set(declared.values())) == 1


@pytest.mark.parametrize("version", ["0.7.0", "v0.7.0"])
def test_normalized_version(version: str) -> None:
    """Release checks accept plain and Git-style version strings."""
    assert normalized_version(version) == "0.7.0"


def test_cleaner_rejects_repository_root() -> None:
    """The cleaning helper cannot recursively remove the repository."""
    with pytest.raises(ValueError, match="outside the repository"):
        safe_directory(ROOT)


def test_cleaner_accepts_scoped_directory() -> None:
    """A generated directory inside the repository is a valid clean target."""
    build_directory = ROOT / "build" / "workflow-test"
    assert safe_directory(build_directory) == Path(build_directory).resolve()
