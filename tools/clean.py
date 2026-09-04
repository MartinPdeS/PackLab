#!/usr/bin/env python3
"""Safely remove PackLab build and manuscript-generated files."""

from __future__ import annotations

import argparse
from pathlib import Path
import shutil


ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT_SUFFIXES = (
    ".aux",
    ".bbl",
    ".blg",
    ".fdb_latexmk",
    ".fls",
    ".glg",
    ".glo",
    ".gls",
    ".ist",
    ".log",
    ".out",
    ".spl",
)


def safe_directory(path: Path) -> Path:
    """Resolve and validate a directory before recursive deletion."""
    resolved = path.resolve()
    if resolved == ROOT or ROOT not in resolved.parents:
        raise ValueError(f"refusing to remove directory outside the repository: {resolved}")
    return resolved


def clean_manuscript(manuscript_dir: Path) -> None:
    """Remove only known LaTeX products from a manuscript directory."""
    directory = safe_directory(manuscript_dir)
    for suffix in MANUSCRIPT_SUFFIXES:
        path = directory / f"packlab{suffix}"
        if path.exists():
            path.unlink()
            print(f"removed {path.relative_to(ROOT)}")
    pdf = directory / "packlab.pdf"
    if pdf.exists():
        pdf.unlink()
        print(f"removed {pdf.relative_to(ROOT)}")


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", type=Path, action="append")
    parser.add_argument(
        "--manuscript-dir",
        type=Path,
        default=ROOT / "docs" / "manuscript" / "softwareX",
    )
    parser.add_argument("--manuscript-only", action="store_true")
    return parser.parse_args()


def main() -> int:
    """Remove requested generated outputs."""
    arguments = parse_arguments()
    if not arguments.manuscript_only:
        build_directories = arguments.build_dir or [ROOT / "build", ROOT / ".skbuild"]
        for requested_directory in build_directories:
            build_directory = safe_directory(requested_directory)
            if build_directory.exists():
                shutil.rmtree(build_directory)
                print(f"removed {build_directory.relative_to(ROOT)}")
    clean_manuscript(arguments.manuscript_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
