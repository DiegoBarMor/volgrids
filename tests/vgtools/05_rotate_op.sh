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
conf_just_stacking="DO_SMIF_STACKING=True DO_SMIF_HBA=False DO_SMIF_HBD=False DO_SMIF_HYDROPHOBIC=False DO_SMIF_HYDROPHILIC=False DO_SMIF_APBS=False SAVE_TRIMMING_MASK=False"

mkdir -p $fop
cp $f_interface/* $fop/
python3 volgrids smiffer prot $fop/prot.pdb --config "$conf_just_stacking"
python3 volgrids smiffer rna  $fop/rna.pdb  --config "$conf_just_stacking"

python3 volgrids vgtools op abs $fop/prot.cmap $fop/prot.abs.cmap
python3 volgrids vgtools op abs $fop/rna.cmap  $fop/rna.abs.cmap
python3 volgrids vgtools op add $fop/prot.cmap $fop/rna.cmap  $fop/prot_plus_rna.cmap
python3 volgrids vgtools op sub $fop/prot.cmap $fop/rna.cmap  $fop/prot_minus_rna.cmap
python3 volgrids vgtools op mul $fop/prot.cmap $fop/rna.cmap  $fop/prot_mul_rna.cmap
python3 volgrids vgtools op div $fop/prot.cmap $fop/rna.cmap  $fop/prot_div_rna.cmap
