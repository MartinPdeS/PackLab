#!/usr/bin/env python3
"""Check the SoftwareX manuscript's reproducibility and submission limits."""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
REQUIRED_SECTIONS = (
    "Motivation and significance",
    "Software description",
    "Illustrative examples",
    "Impact",
    "Conclusions",
)
REQUIRED_DATA = (
    "mctpy_0_0_2_structure_factors.csv",
    "tsang_2001_figure_8_1_3_digitized.csv",
    "tsang_2001_figure_8_3_4_digitized.csv",
)
MAXIMUM_WORDS = 3000
MAXIMUM_FIGURES = 6
MAXIMUM_KEYWORDS = 6
MINIMUM_HIGHLIGHTS = 3
MAXIMUM_HIGHLIGHTS = 5
MAXIMUM_HIGHLIGHT_LENGTH = 85


def texcount(source: Path) -> int:
    """Return the manuscript word count reported by texcount."""
    try:
        completed = subprocess.run(
            ["texcount", "-sum", "-1", source.name],
            cwd=source.parent,
            check=True,
            capture_output=True,
            text=True,
        )
    except FileNotFoundError as error:
        raise RuntimeError("texcount is required for manuscript checking") from error
    except subprocess.CalledProcessError as error:
        raise RuntimeError(f"texcount failed: {error.stderr.strip()}") from error

    try:
        return int(completed.stdout.strip().splitlines()[-1])
    except (IndexError, ValueError) as error:
        raise RuntimeError("could not parse texcount output") from error


def extract_environment(source: str, environment: str) -> str:
    """Return the contents of one LaTeX environment."""
    match = re.search(
        rf"\\begin\{{{re.escape(environment)}\}}(.*?)\\end\{{{re.escape(environment)}\}}",
        source,
        flags=re.DOTALL,
    )
    if match is None:
        raise RuntimeError(f"missing {environment!r} environment")
    return match.group(1)


def inspect_log(log_path: Path) -> list[str]:
    """Return fatal or unresolved-reference diagnostics from a LaTeX log."""
    if not log_path.exists():
        return [f"missing build log: {log_path.relative_to(ROOT)}"]
    log = log_path.read_text(encoding="utf-8", errors="replace")
    patterns = (
        r"Citation .* undefined",
        r"Reference .* undefined",
        r"There were undefined (?:citations|references)",
        r"Undefined control sequence",
        r"Fatal error",
    )
    return [f"LaTeX log contains: {pattern}" for pattern in patterns if re.search(pattern, log)]


def check_manuscript(manuscript_dir: Path) -> tuple[list[str], dict[str, int]]:
    """Return validation failures and summary statistics."""
    source_path = manuscript_dir / "packlab.tex"
    highlights_path = manuscript_dir / "highlights.txt"
    pdf_path = manuscript_dir / "packlab.pdf"
    if not source_path.exists():
        return [f"missing manuscript source: {source_path}"], {}

    source = source_path.read_text(encoding="utf-8")
    failures: list[str] = []
    word_count = texcount(source_path)
    figure_count = len(re.findall(r"\\begin\{figure\}", source))

    try:
        keywords = extract_environment(source, "keyword")
        keyword_count = len(re.split(r"\\sep", keywords))
    except RuntimeError as error:
        failures.append(str(error))
        keyword_count = 0

    highlights = []
    if highlights_path.exists():
        highlights = [
            line.strip() for line in highlights_path.read_text().splitlines() if line.strip()
        ]
    else:
        failures.append("missing highlights.txt")

    if word_count > MAXIMUM_WORDS:
        failures.append(f"word count is {word_count}; limit is {MAXIMUM_WORDS}")
    if figure_count > MAXIMUM_FIGURES:
        failures.append(f"figure count is {figure_count}; limit is {MAXIMUM_FIGURES}")
    if keyword_count > MAXIMUM_KEYWORDS:
        failures.append(f"keyword count is {keyword_count}; limit is {MAXIMUM_KEYWORDS}")
    if highlights and not MINIMUM_HIGHLIGHTS <= len(highlights) <= MAXIMUM_HIGHLIGHTS:
        failures.append(
            f"highlight count is {len(highlights)}; expected "
            f"{MINIMUM_HIGHLIGHTS}--{MAXIMUM_HIGHLIGHTS}"
        )
    for line_number, highlight in enumerate(highlights, start=1):
        if len(highlight) > MAXIMUM_HIGHLIGHT_LENGTH:
            failures.append(
                f"highlight {line_number} has {len(highlight)} characters; "
                f"limit is {MAXIMUM_HIGHLIGHT_LENGTH}"
            )

    for section in REQUIRED_SECTIONS:
        if rf"\section{{{section}}}" not in source:
            failures.append(f"missing required section: {section}")
    identifiers = (
        *[f"C{number}" for number in range(1, 10)],
        *[f"S{number}" for number in range(1, 9)],
    )
    for identifier in identifiers:
        if not re.search(rf"^{identifier}\s*&", source, flags=re.MULTILINE):
            failures.append(f"missing metadata row: {identifier}")

    code_version_match = re.search(
        r"^C1\s*&\s*Current code version\s*&\s*([^\\\s]+)", source, flags=re.MULTILINE
    )
    software_version_match = re.search(
        r"^S1\s*&\s*Current software version\s*&\s*([^\\\s]+)", source, flags=re.MULTILINE
    )
    if code_version_match and software_version_match:
        code_version = code_version_match.group(1)
        software_version = software_version_match.group(1)
        if code_version != software_version:
            failures.append(
                f"code metadata version {code_version} differs from software version "
                f"{software_version}"
            )
        expected_fragments = (f"/tree/v{code_version}", f"/PackLab/{software_version}/")
        for fragment in expected_fragments:
            if fragment not in source:
                failures.append(f"version-specific manuscript link is missing {fragment!r}")
    for filename in REQUIRED_DATA:
        if not (manuscript_dir / "data" / filename).exists():
            failures.append(f"missing reference data: data/{filename}")
    if not pdf_path.exists():
        failures.append("missing compiled manuscript: packlab.pdf")
    failures.extend(inspect_log(manuscript_dir / "packlab.log"))

    statistics = {
        "words": word_count,
        "figures": figure_count,
        "keywords": keyword_count,
        "highlights": len(highlights),
    }
    return failures, statistics


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manuscript-dir",
        type=Path,
        default=ROOT / "docs" / "manuscript" / "softwareX",
    )
    return parser.parse_args()


def main() -> int:
    """Check the manuscript and print an actionable summary."""
    arguments = parse_arguments()
    manuscript_dir = arguments.manuscript_dir.resolve()
    try:
        failures, statistics = check_manuscript(manuscript_dir)
    except RuntimeError as error:
        print(f"manuscript check failed: {error}", file=sys.stderr)
        return 1

    if statistics:
        print(
            "SoftwareX manuscript: "
            f"{statistics['words']} words, {statistics['figures']} figures, "
            f"{statistics['keywords']} keywords, {statistics['highlights']} highlights"
        )
    if failures:
        for failure in failures:
            print(f"ERROR: {failure}", file=sys.stderr)
        return 1

    print("SoftwareX manuscript checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
