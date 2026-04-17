#!/bin/bash
set -eu

echo
echo ">>> TEST VGTOOLS 1: Packing and Unpacking"

folder=testdata/vgtools
fp="$folder/packing"
fu="$folder/unpacking"

############################# PACKING
paths_in="$fp/2esj.hba.mrc $fp/2esj.hbd.mrc $fp/2esj.phi.mrc $fp/2esj.pho.mrc $fp/2esj.stk.mrc"
path_out="$fp/2esj.cmap"

# shellcheck disable=SC2086
python3 volgrids vgtools pack -i $paths_in -o "$path_out"


############################# UNPACKING
path_in="$fu/1iqj.cmap"
python3 volgrids vgtools unpack -i "$path_in"
