#!/usr/bin/env python
"""Build the documentation for FREQCOR.

Run this script from the project root directory.
"""

import subprocess
import sys
import webbrowser
import os
from pathlib import Path


def _rmtree_onerror(func, path, exc_info):
    """Retry removing read-only files/dirs on Windows.

    Sphinx sometimes leaves read-only files in the build directory; on Windows this
    can cause shutil.rmtree to fail with PermissionError.
    """
    try:
        import os
        import stat

        os.chmod(path, stat.S_IWRITE)
        func(path)
    except Exception:
        raise


def build_docs():
    """Build the Sphinx documentation."""
    docs_dir = Path(__file__).parent / 'docs'
    build_dir = docs_dir / 'build' / 'html'

    if build_dir.exists():
        import shutil
        shutil.rmtree(build_dir, onerror=_rmtree_onerror)

    build_dir.parent.mkdir(parents=True, exist_ok=True)

    source_dir = docs_dir / 'source'
    build_root = docs_dir / 'build'

    # Use the current Python interpreter to run Sphinx. This avoids relying on a
    # sphinx-build executable being on PATH (common issue on Windows).
    result = subprocess.run(
        [
            sys.executable,
            '-m',
            'sphinx',
            '-M',
            'html',
            str(source_dir),
            str(build_root),
        ],
        cwd=docs_dir,
    )

    if result.returncode != 0:
        return False

    index_path = build_dir / 'index.html'
    if index_path.exists():
        if not (os.environ.get('CI') or os.environ.get('GITHUB_ACTIONS')):
            webbrowser.open(index_path.as_uri())
        return True

    return False


if __name__ == '__main__':
    raise SystemExit(0 if build_docs() else 1)
