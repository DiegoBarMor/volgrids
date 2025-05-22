#!/bin/bash
set -eu

pip install .
rm -rf build volpot.egg-info
