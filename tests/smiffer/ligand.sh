#!/bin/bash
set -eu

echo
echo ">>> TEST 4: Ligands"

folder=testdata/smiffer/ligand

### Ligand: Phenol
python3 smiffer.py ligand $folder/1lph.pdb -b $folder/phenol.chem

mv "$folder/1lph.cmap" "$folder/phenol.cmap"
mv "$folder/1lph.meta.json" "$folder/phenol.meta.json"

### Protein system (for reference)
python3 smiffer.py prot $folder/1lph.pdb
