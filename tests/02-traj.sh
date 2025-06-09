#!/bin/bash
set -eu

echo
echo ">>> TEST 02: Trajectory mode"

folder=data/tests/02-traj
python3 smiffer.py rna $folder/7vki.pdb -o $folder -t $folder/7vki.xtc
rm -f $folder/.[!.]*.npz $folder/.[!.]*.lock
