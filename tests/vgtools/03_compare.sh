#!/bin/bash
set -eu

echo
echo ">>> TEST VGTOOLS 3: Comparison of grids"

folder="testdata/vgtools"
fc="$folder/converting"
fp="$folder/packing"

############################# COMPARE BETWEEN FORMATS
path_dx="$fc/1iqj.stk.dx"
path_mrc="$fc/1iqj.stk.mrc"
path_ccp4="$fc/1iqj.stk.ccp4"
path_cmap="$fc/1iqj.stk.cmap"

python3 volgrids vgtools compare "$path_dx" "$path_dx"   -t 1e-4; echo
python3 volgrids vgtools compare "$path_dx" "$path_mrc"  -t 1e-4; echo
python3 volgrids vgtools compare "$path_dx" "$path_ccp4" -t 1e-4; echo
python3 volgrids vgtools compare "$path_dx" "$path_cmap" -t 1e-4; echo

python3 volgrids vgtools compare "$path_mrc" "$path_mrc"  -t 1e-20; echo
python3 volgrids vgtools compare "$path_mrc" "$path_ccp4" -t 1e-20; echo
python3 volgrids vgtools compare "$path_mrc" "$path_cmap" -t 1e-20; echo

python3 volgrids vgtools compare "$path_ccp4" "$path_ccp4" -t 1e-20; echo
python3 volgrids vgtools compare "$path_ccp4" "$path_cmap" -t 1e-20; echo

python3 volgrids vgtools compare "$path_cmap" "$path_cmap" -t 1e-20; echo


############################# DIFFERENT THRESHOLDS
python3 volgrids vgtools compare "$path_dx" "$path_mrc" -t 0     ; echo
python3 volgrids vgtools compare "$path_dx" "$path_mrc" -t 1e-20 ; echo
python3 volgrids vgtools compare "$path_dx" "$path_mrc" -t 1e-10 ; echo
python3 volgrids vgtools compare "$path_dx" "$path_mrc" -t 1e-5  ; echo
python3 volgrids vgtools compare "$path_dx" "$path_mrc" -t 1e-4  ; echo
python3 volgrids vgtools compare "$path_dx" "$path_mrc" -t 1e-3  ; echo
python3 volgrids vgtools compare "$path_dx" "$path_mrc" -t 1e-2  ; echo
python3 volgrids vgtools compare "$path_dx" "$path_mrc" -t 1e-1  ; echo
python3 volgrids vgtools compare "$path_dx" "$path_mrc" -t 1     ; echo


############################# INCOMPATIBLE GRIDS (DIFFERENT RESOLUTIONS)
path_grid0="$fc/1iqj.stk.mrc"
path_grid1="$fp/2esj.stk.mrc"

python3 volgrids vgtools compare "$path_grid0" "$path_grid1"
