#!/bin/bash
set -eu

pip uninstall -y volgrids
pip install .
rm -rf build volgrids.egg-info
