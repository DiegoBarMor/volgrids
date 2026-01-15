#!/bin/bash
set -eu

echo
echo ">>> TEST SMIFFER 3: Trajectory mode"

folder_trajs=testdata/smiffer/trajs
folder_few_frames=$folder_trajs/few_frames
folder_few_resids=$folder_trajs/few_resids

python3 volgrids smiffer rna $folder_few_resids/7vki.pdb -o $folder_few_resids -t $folder_few_resids/7vki.xtc
python3 volgrids smiffer rna $folder_few_frames/fse.pdb  -o $folder_few_frames -t $folder_few_frames/fse.xtc -s 43.66 38.12 36.13 12

rm -f $folder_few_frames/.[!.]*.npz $folder_few_frames/.[!.]*.lock
rm -f $folder_few_resids/.[!.]*.npz $folder_few_resids/.[!.]*.lock
