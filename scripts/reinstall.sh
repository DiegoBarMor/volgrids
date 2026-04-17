#!/bin/bash
set -eu

### Re-install the package locally

pip uninstall volgrids -y || true
pip install .
rm -rf build volgrids.egg-info
