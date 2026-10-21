"""Sphinx configuration for PackLab documentation."""

import os
import sys
from importlib.metadata import version as distribution_version
from pathlib import Path

import matplotlib.pyplot as plt

from martinpdes_docs import configure


package_name = "PackLab"
version = os.getenv("tag", distribution_version(package_name))

# The documentation workflow installs the built wheel before invoking Sphinx.
# Keep the gallery examples on that wheel so they use PackLab's native core.
sys.path.insert(0, str(Path(__file__).resolve().parent))
os.environ["PACKLAB_DOCS_ROOT"] = str(Path(__file__).resolve().parents[1])


def reset_mpl(gallery_conf, fname):
    """Restore Matplotlib defaults between gallery examples."""
    plt.rcdefaults()


globals().update(
    configure(
        project=package_name,
        package=package_name,
        repository="MartinPdeS/PackLab",
        logo="_static/packlab-mark.png",
        favicon="_static/packlab-favicon.png",
        prefer_installed_package=True,
        gallery_default_thumb="_static/gallery-default.png",
        gallery_options={
            "reset_modules": "conf.reset_mpl",
        },
        html_theme_options={
            "logo": {"text": package_name, "image": "_static/packlab-mark.png"},
            "show_nav_level": 1,
            "navbar_end": ["version-switcher", "theme-switcher", "navbar-icon-links"],
            "show_prev_next": True,
        },
    )
)

# PackLab keeps its responsive article and navigation refinements in these
# project-owned stylesheets; the shared stylesheet remains loaded as well.
html_css_files = ["default.css", "layout.css", "martinpdes_layout.css"]

autodoc_default_options.update(
    {
        "members": True,
        "member-order": "bysource",
    }
)
html_context.update({"doc_path": "docs/source", "default_mode": "dark"})

major, minor = version.split(".")[:2]
binder_branch = f"v{major}.{minor}.x"
