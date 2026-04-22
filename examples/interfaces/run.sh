#!/bin/bash
set -euo pipefail

fdata="testdata/smiffer/interfaces/prot_rna"
conf_just_stacking="DO_SMIF_STACKING=True DO_SMIF_HBA=False DO_SMIF_HBD=False DO_SMIF_HYDROPHOBIC=False DO_SMIF_HYDROPHILIC=False DO_SMIF_APBS=False SAVE_TRIMMING_MASK=False"

python3 volgrids smiffer prot $fdata/prot.pdb -c "$conf_just_stacking"
python3 volgrids smiffer rna  $fdata/rna.pdb  -c "$conf_just_stacking"
mv $fdata/prot.cmap $fdata/prot.smif.cmap
mv $fdata/rna.cmap $fdata/rna.smif.cmap

python3 volgrids smutils occupancy prot $fdata/prot.pdb -c "$conf_just_stacking"
python3 volgrids smutils occupancy rna  $fdata/rna.pdb -c "$conf_just_stacking"
mv $fdata/prot.cmap $fdata/prot.og.cmap
mv $fdata/rna.cmap $fdata/rna.og.cmap

# [WIP]
