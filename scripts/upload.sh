#!/bin/bash
set -eu

### Build and upload the package to PyPI

# pip install build twine
python3 -m build
twine upload dist/*
rm -rf build dist volgrids.egg-info
