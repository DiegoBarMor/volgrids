#!/bin/bash
set -eu

echo
echo ">>> TEST SMIFFER 3: Trajectory mode"

folder=testdata/smiffer/traj
python3 run/smiffer.py rna $folder/7vki.pdb -o $folder -t $folder/7vki.xtc
rm -f $folder/.[!.]*.npz $folder/.[!.]*.lock
