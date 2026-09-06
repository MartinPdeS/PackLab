#!/usr/bin/env python3
"""Print the next semantic-version tag from the latest reachable release tag."""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
TAG_PATTERN = re.compile(r"v(?P<major>0|[1-9]\d*)\.(?P<minor>0|[1-9]\d*)\.(?P<patch>0|[1-9]\d*)$")


def latest_version() -> tuple[int, int, int]:
    """Return the latest reachable semantic-version tag."""
    completed = subprocess.run(
        ["git", "describe", "--tags", "--abbrev=0", "--match", "v[0-9]*"],
        cwd=ROOT,
        check=False,
        text=True,
        capture_output=True,
    )
    if completed.returncode != 0:
        raise RuntimeError("no reachable semantic-version tag was found")
    match = TAG_PATTERN.fullmatch(completed.stdout.strip())
    if match is None:
        raise RuntimeError("latest reachable tag is not vMAJOR.MINOR.PATCH")
    return tuple(int(match.group(name)) for name in ("major", "minor", "patch"))


def main() -> int:
    """Print the requested semantic-version increment."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("kind", choices=("major", "minor", "patch"))
    arguments = parser.parse_args()
    try:
        major, minor, patch = latest_version()
    except RuntimeError as error:
        print(f"could not determine next release: {error}", file=sys.stderr)
        return 1

    if arguments.kind == "major":
        major, minor, patch = major + 1, 0, 0
    elif arguments.kind == "minor":
        minor, patch = minor + 1, 0
    else:
        patch += 1
    print(f"v{major}.{minor}.{patch}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
