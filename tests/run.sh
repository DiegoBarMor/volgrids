#!/bin/bash
set -eu

bash tests/00-pocket_sphere.sh
bash tests/01-whole.sh
python3 tests/02-conversions.py
python3 tests/03-cmap.py
