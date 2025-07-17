#!/bin/bash
set -eu

echo
echo ">>> TEST SMIFFER 4: Ligands"

folder=testdata/smiffer/ligand

### Ligand: Phenol
python3 run/smiffer.py ligand $folder/1lph.pdb -b $folder/phenol.chem
mv "$folder/1lph.cmap" "$folder/phenol.cmap"

### Protein system (for reference)
python3 run/smiffer.py prot $folder/1lph.pdb
