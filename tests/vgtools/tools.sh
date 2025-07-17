#!/bin/bash
set -eu

echo
echo ">>> TEST VGTOOLS 0: Conversions between DX,MRC,CCP4,CMAP formats; packing; unpacking"

folder=testdata/vgtools
fc="$folder/converting"
fp="$folder/packing"
fu="$folder/unpacking"
ff="$folder/fix_cmap"

############################# CONVERSIONS
path_dx_input="$fc/1iqj.stk.dx"
path_mrc_input="$fc/1iqj.stk.mrc"
path_ccp4_input="$fc/1iqj.stk.ccp4"
path_cmap_input="$fc/1iqj.stk.cmap"

path_dx_to_dx="$fc/dx-dx.dx"
path_dx_to_mrc="$fc/dx-mrc.mrc"
path_dx_to_ccp4="$fc/dx-ccp4.ccp4"
path_dx_to_cmap="$fc/dx-cmap.cmap"

path_mrc_to_dx="$fc/mrc-dx.dx"
path_mrc_to_mrc="$fc/mrc-mrc.mrc"
path_mrc_to_ccp4="$fc/mrc-ccp4.ccp4"
path_mrc_to_cmap="$fc/mrc-cmap.cmap"

path_ccp4_to_dx="$fc/ccp4-dx.dx"
path_ccp4_to_mrc="$fc/ccp4-mrc.mrc"
path_ccp4_to_ccp4="$fc/ccp4-ccp4.ccp4"
path_ccp4_to_cmap="$fc/ccp4-cmap.cmap"

path_cmap_to_dx="$fc/cmap-dx.dx"
path_cmap_to_mrc="$fc/cmap-mrc.mrc"
path_cmap_to_ccp4="$fc/cmap-ccp4.ccp4"
path_cmap_to_cmap="$fc/cmap-cmap.cmap"

python3 run/vgtools.py convert "$path_dx_input"   \
    --dx "$path_dx_to_dx"       --mrc "$path_dx_to_mrc"   \
    --ccp4 "$path_dx_to_ccp4"   --cmap "$path_dx_to_cmap"

python3 run/vgtools.py convert "$path_mrc_input"  \
    --dx "$path_mrc_to_dx"      --mrc "$path_mrc_to_mrc"  \
    --ccp4 "$path_mrc_to_ccp4"  --cmap "$path_mrc_to_cmap"

python3 run/vgtools.py convert "$path_ccp4_input" \
    --dx "$path_ccp4_to_dx"     --mrc "$path_ccp4_to_mrc" \
    --ccp4 "$path_ccp4_to_ccp4" --cmap "$path_ccp4_to_cmap"

python3 run/vgtools.py convert "$path_cmap_input" \
    --dx "$path_cmap_to_dx"     --mrc "$path_cmap_to_mrc" \
    --ccp4 "$path_cmap_to_ccp4" --cmap "$path_cmap_to_cmap"



############################# PACKING
paths_in="$fp/2esj.hba.mrc $fp/2esj.hbd.mrc $fp/2esj.phi.mrc $fp/2esj.pho.mrc $fp/2esj.stk.mrc"
path_out="$fp/2esj.cmap"

# shellcheck disable=SC2086
python3 run/vgtools.py pack -i $paths_in -o "$path_out"


############################# UNPACKING
path_in="$fu/1iqj.cmap"
python3 run/vgtools.py unpack -i "$path_in"


############################# FIX CMAP
path_in="$ff/hbdonors.issue.cmap"
path_out="$ff/hbdonors.fixed.cmap"
python3 run/vgtools.py fix_cmap -i "$path_in" -o "$path_out"
