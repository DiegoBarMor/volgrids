#!/bin/bash
set -eu

echo
echo ">>> TEST VGTOOLS 0: Conversions between MRC,DX,CMAP formats; packing; unpacking"

folder=testdata/vgtools
fc="$folder/converting"
fp="$folder/packing"
fu="$folder/unpacking"
ff="$folder/fix_cmap"

############################# CONVERSIONS
path_dx_input="$fc/1iqj.stk.dx"
path_mrc_input="$fc/1iqj.stk.mrc"
path_cmap_input="$fc/1iqj.stk.cmap"

path_dx_to_dx="$fc/dx2dx.dx"
path_dx_to_mrc="$fc/dx2mrc.mrc"
path_dx_to_cmap="$fc/dx2cmap.cmap"

path_mrc_to_dx="$fc/mrc2dx.dx"
path_mrc_to_mrc="$fc/mrc2mrc.mrc"
path_mrc_to_cmap="$fc/mrc2cmap.cmap"

path_cmap_to_dx="$fc/cmap2dx.dx"
path_cmap_to_mrc="$fc/cmap2mrc.mrc"
path_cmap_to_cmap="$fc/cmap2cmap.cmap"

python3 vgtools.py convert "$path_dx_input"   --dx "$path_dx_to_dx"   --mrc "$path_dx_to_mrc"   --cmap "$path_dx_to_cmap"
python3 vgtools.py convert "$path_mrc_input"  --dx "$path_mrc_to_dx"  --mrc "$path_mrc_to_mrc"  --cmap "$path_mrc_to_cmap"
python3 vgtools.py convert "$path_cmap_input" --dx "$path_cmap_to_dx" --mrc "$path_cmap_to_mrc" --cmap "$path_cmap_to_cmap"


############################# PACKING
paths_in="$fp/2esj.hba.mrc $fp/2esj.hbd.mrc $fp/2esj.phi.mrc $fp/2esj.pho.mrc $fp/2esj.stk.mrc"
path_out="$fp/2esj.cmap"

# shellcheck disable=SC2086
python3 vgtools.py pack -i $paths_in -o "$path_out"


############################# UNPACKING
path_in="$fu/1iqj.cmap"
python3 vgtools.py unpack -i "$path_in"


############################# FIX CMAP
path_in="$ff/hbdonors.issue.cmap"
path_out="$ff/hbdonors.fixed.cmap"
python3 vgtools.py fix-cmap -i "$path_in" -o "$path_out"
