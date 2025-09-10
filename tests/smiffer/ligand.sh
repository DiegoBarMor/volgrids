#!/bin/bash
set -eu

echo
echo ">>> TEST SMIFFER 4: Ligands"

folder=testdata/smiffer/ligands

### Ligand: TSC aka (1S)-1-AMINO-2-(1H-INDOL-3-YL)ETHANOL
python3 run/smiffer.py ligand $folder/all_interactions/tsc.pdb     -b $folder/all_interactions/tsc.chem
python3 run/smiffer.py ligand $folder/all_interactions/tsc_noH.pdb -b $folder/all_interactions/tsc.chem

### Ligand: 9U3 aka N,N,N',N'-tetramethylethane-1,2-diamine
python3 run/smiffer.py ligand $folder/nr3_case/9u3.pdb -b $folder/nr3_case/9u3.chem

# ### Protein system (for reference)
# python3 run/smiffer.py prot $folder/1lph.pdb
