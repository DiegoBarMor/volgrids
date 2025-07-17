#!/bin/bash
set -eu

pip uninstall -y volgrids
pip install .
rm -rf build src/volgrids.egg-info
