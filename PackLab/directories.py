#!/usr/bin/env python
# -*- coding: utf-8; py-indent-offset:4 -*-

from pathlib import Path
import PackLab

root_path = Path(PackLab.__path__[0])

repository = project_path = Path(root_path)

doc_css_path = project_path.parent / "docs/source/_static/default.css"

data = repository / "data"

