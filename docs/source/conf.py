"""Sphinx configuration for PackLab documentation."""

import os
from pathlib import Path

from martinpdes_docs import configure


package_name = "PackLab"
docs_source = Path(__file__).resolve().parent

# Validation examples use this root to locate manuscript reference data.
os.environ["PACKLAB_DOCS_ROOT"] = str(docs_source.parent)

globals().update(
    configure(
        project=package_name,
        package=package_name,
        version="auto",
        author="Martin Poinsinet de Sivry-Houle",
        repository=f"MartinPdeS/{package_name}",
        logo="../images/logo.svg",
        favicon="_static/packlab-favicon.png",
        prefer_installed_package=True,
        docs_root=docs_source.parent,
        profile="scientific-gallery",
        gallery_default_thumb=str(docs_source / "_static" / "gallery-default.png"),
        gallery_options={"capture_repr": ()},
        logo_text="",
        logo_alt_text=f"{package_name} hard-sphere packing logo",
        html_css_files=["layout.css"],
        autodoc_default_options={
            "members": True,
            "member-order": "bysource",
        },
    )
)
