#!/bin/bash
set -eu

echo
echo ">>> TEST SMIFFER 4: Ligands"

folder=testdata/smiffer/ligands

### Ligand: TSC [(1S)-1-AMINO-2-(1H-INDOL-3-YL)ETHANOL]
python3 run/smiffer.py ligand $folder/tsc/tsc.pdb     -b $folder/tsc/params/tsc.chem
python3 run/smiffer.py ligand $folder/tsc/tsc_noH.pdb -b $folder/tsc/params/tsc.chem


# ### Protein system (for reference)
# python3 run/smiffer.py prot $folder/1lph.pdb
