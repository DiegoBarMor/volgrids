#!/bin/bash
set -euo pipefail

fdata="testdata/smiffer/interfaces/prot_rna"
config="examples/interfaces/config.ini"

python3 volgrids smiffer prot $fdata/prot.pdb -c "$config"
python3 volgrids smiffer rna  $fdata/rna.pdb  -c "$config"
mv $fdata/prot.cmap $fdata/prot.smif.cmap
mv $fdata/rna.cmap $fdata/rna.smif.cmap

python3 volgrids smutils occupancy prot $fdata/prot.pdb -c "$config"
python3 volgrids smutils occupancy rna  $fdata/rna.pdb  -c "$config"
mv $fdata/prot.cmap $fdata/prot.og.cmap
mv $fdata/rna.cmap $fdata/rna.og.cmap

python3 volgrids vgtools op mul $fdata/prot.smif.cmap $fdata/rna.og.cmap $fdata/prot_smif.rna_og.cmap
python3 volgrids vgtools op mul $fdata/rna.smif.cmap $fdata/prot.og.cmap $fdata/rna_smif.prot_og.cmap
