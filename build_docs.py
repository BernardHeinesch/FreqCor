#!/usr/bin/env python
"""Build the documentation for FREQCOR.

Run this script from the project root directory.
"""

import subprocess
import sys
import webbrowser
from pathlib import Path


def build_docs():
    """Build the Sphinx documentation."""
    docs_dir = Path(__file__).parent / 'docs'
    build_dir = docs_dir / 'build' / 'html'

    if build_dir.exists():
        import shutil
        shutil.rmtree(build_dir)

    build_dir.parent.mkdir(parents=True, exist_ok=True)

    if sys.platform.startswith('win'):
        result = subprocess.run(['cmd', '/c', 'make.bat', 'html'], cwd=docs_dir)
    else:
        result = subprocess.run(['make', 'html'], cwd=docs_dir)

    if result.returncode != 0:
        return False

    index_path = build_dir / 'index.html'
    if index_path.exists():
        webbrowser.open(index_path.as_uri())
        return True

    return False


if __name__ == '__main__':
    raise SystemExit(0 if build_docs() else 1)
