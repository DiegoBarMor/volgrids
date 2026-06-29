#!/bin/bash
set -euo pipefail

echo
echo ">>> TEST VGTOOLS 5: Rotation and operations (abs, add, sub, mul, div) on grids"

f_interface="testdata/smiffer/interfaces/prot_rna"

folder="testdata/vgtools"
# frot="$folder/rotate"
fop="$folder/operations"

############################# GRID ROTATIONS
# [TODO]


############################# OPERATIONS BETWEEN GRIDS
conf_just_stacking="SMIF_STK=True SMIF_HBA=False SMIF_HBD=False SMIF_HPHOB=False SMIF_HPHIL=False SMIF_APBS=False TRIM_SAVE=False"

mkdir -p $fop
cp $f_interface/* $fop/
python3 volgrids smiffer $fop/prot.pdb --config "$conf_just_stacking"
python3 volgrids smiffer $fop/rna.pdb  --config "$conf_just_stacking"

python3 volgrids vgtools op abs $fop/prot.stk.mrc $fop/prot.abs.cmap
python3 volgrids vgtools op abs $fop/rna.stk.mrc  $fop/rna.abs.cmap
python3 volgrids vgtools op add $fop/prot.stk.mrc $fop/rna.stk.mrc  $fop/prot_plus_rna.cmap
python3 volgrids vgtools op sub $fop/prot.stk.mrc $fop/rna.stk.mrc  $fop/prot_minus_rna.cmap
python3 volgrids vgtools op mul $fop/prot.stk.mrc $fop/rna.stk.mrc  $fop/prot_mul_rna.cmap
python3 volgrids vgtools op div $fop/prot.stk.mrc $fop/rna.stk.mrc  $fop/prot_div_rna.cmap
