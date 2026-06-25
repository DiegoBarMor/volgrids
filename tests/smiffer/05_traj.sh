#!/bin/bash
set -eu

echo
echo ">>> TEST SMIFFER 5: Trajectory mode"

folder_trajs=testdata/smiffer/trajs
folder_few_resids=$folder_trajs/few_resids

python3 volgrids smiffer $folder_few_resids/7vki.pdb -o $folder_few_resids -t $folder_few_resids/7vki.xtc
