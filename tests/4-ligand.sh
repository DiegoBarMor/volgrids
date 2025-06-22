#!/bin/bash
set -eu

echo
echo ">>> TEST 4: Ligands"

folder=data/tests/4-ligand

### Ligand: Phenol
python3 smiffer.py ligand $folder/1lph.pdb -b $folder/phenol.atoms

mapfile -t paths_grids < <(ls "$folder/1lph".*.cmap) # get the list of output files
python3 vgtools.py pack -i "${paths_grids[@]}" -o "$folder/phenol.cmap"  # pack them into a single CMAP file
rm -f "${paths_grids[@]}" # remove the individual grid files
mv "$folder/1lph.meta.json" "$folder/phenol.meta.json"

### Protein system (for reference)
python3 smiffer.py prot $folder/1lph.pdb

mapfile -t paths_grids < <(ls "$folder/1lph".*.cmap)
python3 vgtools.py pack -i "${paths_grids[@]}" -o "$folder/1lph.cmap"
rm -f "${paths_grids[@]}"
