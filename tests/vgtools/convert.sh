#!/bin/bash
set -eu

echo
echo ">>> TEST VGTOOLS 0: Conversions between DX,MRC,CCP4,CMAP formats"

folder="testdata/vgtools/converting"

path_dx_input="$folder/1iqj.stk.dx"
path_mrc_input="$folder/1iqj.stk.mrc"
path_ccp4_input="$folder/1iqj.stk.ccp4"
path_cmap_input="$folder/1iqj.stk.cmap"

path_dx_to_dx="$folder/dx-dx.dx"
path_dx_to_mrc="$folder/dx-mrc.mrc"
path_dx_to_ccp4="$folder/dx-ccp4.ccp4"
path_dx_to_cmap="$folder/dx-cmap.cmap"

path_mrc_to_dx="$folder/mrc-dx.dx"
path_mrc_to_mrc="$folder/mrc-mrc.mrc"
path_mrc_to_ccp4="$folder/mrc-ccp4.ccp4"
path_mrc_to_cmap="$folder/mrc-cmap.cmap"

path_ccp4_to_dx="$folder/ccp4-dx.dx"
path_ccp4_to_mrc="$folder/ccp4-mrc.mrc"
path_ccp4_to_ccp4="$folder/ccp4-ccp4.ccp4"
path_ccp4_to_cmap="$folder/ccp4-cmap.cmap"

path_cmap_to_dx="$folder/cmap-dx.dx"
path_cmap_to_mrc="$folder/cmap-mrc.mrc"
path_cmap_to_ccp4="$folder/cmap-ccp4.ccp4"
path_cmap_to_cmap="$folder/cmap-cmap.cmap"

python3 vgtools.py convert "$path_dx_input"   \
    --dx "$path_dx_to_dx"       --mrc "$path_dx_to_mrc"   \
    --ccp4 "$path_dx_to_ccp4"   --cmap "$path_dx_to_cmap"

python3 vgtools.py convert "$path_mrc_input"  \
    --dx "$path_mrc_to_dx"      --mrc "$path_mrc_to_mrc"  \
    --ccp4 "$path_mrc_to_ccp4"  --cmap "$path_mrc_to_cmap"

python3 vgtools.py convert "$path_ccp4_input" \
    --dx "$path_ccp4_to_dx"     --mrc "$path_ccp4_to_mrc" \
    --ccp4 "$path_ccp4_to_ccp4" --cmap "$path_ccp4_to_cmap"

python3 vgtools.py convert "$path_cmap_input" \
    --dx "$path_cmap_to_dx"     --mrc "$path_cmap_to_mrc" \
    --ccp4 "$path_cmap_to_ccp4" --cmap "$path_cmap_to_cmap"
