#!/bin/bash
set -eu

echo
echo ">>> TEST SMIFFER 2: Trajectory mode"

folder=testdata/smiffer/traj
python3 smiffer.py rna $folder/7vki.pdb -o $folder -t $folder/7vki.xtc
rm -f $folder/.[!.]*.npz $folder/.[!.]*.lock
