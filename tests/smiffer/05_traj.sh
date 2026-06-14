#!/bin/bash
set -eu

echo
echo ">>> TEST SMIFFER 5: Trajectory mode"

folder_trajs=testdata/smiffer/trajs
folder_few_frames=$folder_trajs/few_frames
folder_few_resids=$folder_trajs/few_resids

python3 volgrids smiffer $folder_few_resids/7vki.pdb -o $folder_few_resids -t $folder_few_resids/7vki.xtc

path_pdb=$folder_few_frames/fse.pdb
path_traj=$folder_few_frames/fse.xtc
query="resid 35 36 37 38"
spheres=$(python3 volgrids smutils sphere find $path_pdb "$query" -t $path_traj -r 2.0)

python3 volgrids smiffer $folder_few_frames/fse.pdb -o $folder_few_frames -t $folder_few_frames/fse.xtc -s "$spheres" -n 4
