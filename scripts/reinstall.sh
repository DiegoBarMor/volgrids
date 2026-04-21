#!/bin/bash
set -eu

### Re-install the package locally

bash scripts/_prepare.sh

pip uninstall volgrids -y || true
pip install .
rm -rf build volgrids.egg-info
