#!/usr/bin/env python3
"""Diagnose the tools and Python packages used by PackLab development."""

from __future__ import annotations

import importlib
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
COMMANDS = ("git", "cmake", "latexmk", "texcount")
MODULES = ("PackLab", "numpy", "pybind11", "pytest", "ruff", "sphinx")
CACHE_ROOT = Path(tempfile.gettempdir()) / "packlab-doctor-cache"
CACHE_ROOT.mkdir(parents=True, exist_ok=True)
MATPLOTLIB_CACHE = CACHE_ROOT / "matplotlib"
MATPLOTLIB_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MATPLOTLIB_CACHE))
os.environ.setdefault("XDG_CACHE_HOME", str(CACHE_ROOT))


def command_version(command: str) -> str | None:
    """Return a concise version line for an executable, if installed."""
    executable = shutil.which(command)
    if executable is None:
        return None
    completed = subprocess.run(
        [executable, "--version"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    output = completed.stdout.strip() or completed.stderr.strip()
    return output.splitlines()[0] if output else executable


def module_version(module_name: str) -> str | None:
    """Return an importable module's reported version."""
    try:
        module = importlib.import_module(module_name)
    except Exception:
        return None
    return str(getattr(module, "__version__", "installed"))


def main() -> int:
    """Print diagnostics and fail when the development environment is incomplete."""
    failures: list[str] = []
    python_ok = sys.version_info >= (3, 10)
    print(f"{'OK' if python_ok else 'MISSING'} Python: {sys.version.split()[0]} ({sys.executable})")
    if not python_ok:
        failures.append("Python 3.10 or newer is required")

    for command in COMMANDS:
        version = command_version(command)
        print(f"{'OK' if version else 'MISSING'} command {command}: {version or 'not found'}")
        if version is None:
            failures.append(f"command not found: {command}")

    compiler = next(
        ((name, command_version(name)) for name in ("c++", "clang++", "g++") if shutil.which(name)),
        None,
    )
    print(
        f"{'OK' if compiler else 'MISSING'} C++ compiler: "
        f"{compiler[1] if compiler else 'not found'}"
    )
    if compiler is None:
        failures.append("a C++ compiler is required")

    for module_name in MODULES:
        version = module_version(module_name)
        print(
            f"{'OK' if version else 'MISSING'} Python module {module_name}: "
            f"{version or 'not importable'}"
        )
        if version is None:
            failures.append(f"Python module not importable: {module_name}")

    if failures:
        print("\nEnvironment needs attention:", file=sys.stderr)
        for failure in failures:
            print(f"- {failure}", file=sys.stderr)
        print("Run `make setup` after installing the missing system tools.", file=sys.stderr)
        return 1

    print("\nPackLab development environment is ready.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
