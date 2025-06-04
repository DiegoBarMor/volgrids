#!/bin/bash
set -eu

echo
echo ">>> TEST 04: Trajectory mode"

folder=data/tests/04-traj
python3 smiffer.py -i $folder/7vki.pdb -o $folder -t $folder/7vki.xtc -n -w
rm -f $folder/.[!.]*.npz $folder/.[!.]*.lock
