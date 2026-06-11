#!/bin/bash
set -eu

echo
echo ">>> TEST VGTOOLS 0: Conversions between DX,MRC,CCP4,CMAP formats"

fdata="testdata/vgtools/converting"
fout_dx="$fdata/out_dx"
fout_mrc="$fdata/out_mrc"
fout_ccp4="$fdata/out_ccp4"
fout_cmap="$fdata/out_cmap"

path_dx_input="$fdata/1iqj.stk.dx"
path_mrc_input="$fdata/1iqj.stk.mrc"
path_ccp4_input="$fdata/1iqj.stk.ccp4"
path_cmap_input="$fdata/1iqj.stk.cmap"

python3 volgrids vgtools convert "$path_dx_input"   "$fout_dx" -f DX
python3 volgrids vgtools convert "$path_dx_input"   "$fout_dx" -f MRC
python3 volgrids vgtools convert "$path_dx_input"   "$fout_dx" -f CCP4
python3 volgrids vgtools convert "$path_dx_input"   "$fout_dx" -f CMAP

python3 volgrids vgtools convert "$path_mrc_input"  "$fout_mrc" -f DX
python3 volgrids vgtools convert "$path_mrc_input"  "$fout_mrc" -f MRC
python3 volgrids vgtools convert "$path_mrc_input"  "$fout_mrc" -f CCP4
python3 volgrids vgtools convert "$path_mrc_input"  "$fout_mrc" -f CMAP

python3 volgrids vgtools convert "$path_ccp4_input" "$fout_ccp4" -f DX
python3 volgrids vgtools convert "$path_ccp4_input" "$fout_ccp4" -f MRC
python3 volgrids vgtools convert "$path_ccp4_input" "$fout_ccp4" -f CCP4
python3 volgrids vgtools convert "$path_ccp4_input" "$fout_ccp4" -f CMAP

python3 volgrids vgtools convert "$path_cmap_input" "$fout_cmap" -f DX
python3 volgrids vgtools convert "$path_cmap_input" "$fout_cmap" -f MRC
python3 volgrids vgtools convert "$path_cmap_input" "$fout_cmap" -f CCP4
python3 volgrids vgtools convert "$path_cmap_input" "$fout_cmap" -f CMAP
