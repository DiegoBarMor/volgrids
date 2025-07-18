#!/bin/bash
set -eu

echo
echo ">>> TEST VGTOOLS 2: Fixing CMAP files with different sized grids"

folder="testdata/vgtools/fix_cmap"

path_in="$folder/hbdonors.issue.cmap"
path_out="$folder/hbdonors.fixed.cmap"
python3 run/vgtools.py fix_cmap -i "$path_in" -o "$path_out"
